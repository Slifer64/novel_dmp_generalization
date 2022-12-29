%% FileIO class
% 

classdef FileIO < matlab.mixin.Copyable

    methods (Access = public)
        
        %% FileIO constructor.
        %  @param[in] filename: Name of the file to open.
        %  @param[in] open_mode: Determines the open mode. Combination of flags from @FileIO.OpenMode. (optional, default=in|out)
        function this = FileIO(filename, open_mode)
            
            if (nargin < 2), open_mode = bitor(FileIO.in, FileIO.out); end
            
            this.setMatType(FileIO.ARMA); % set default matrix type to ARMA
            
            this.header_start = 0;
                        
            this.in_flag = bitand(open_mode, FileIO.in);
            this.out_flag = bitand(open_mode, FileIO.out);
            trunc_flag = bitand(open_mode, FileIO.trunc);

            % check if file exists
            file_exist = exist(filename, 'file') == 2; % replace by isfile in R2017b and above

            if (~this.out_flag)
                if (~this.in_flag), error('[FileIO.FileIO]: Open-mode "in" and/or "out" must be specified!'); end
                if (this.in_flag && trunc_flag), error('[FileIO.FileIO]: Setting only "in" and "trunc" makes no sense :P ...'); end
                if (this.in_flag && ~file_exist), error(['[FileIO.FileIO]: Open-mode "in" was specified but file "' filename '" does not exist or cannot be accessed...']); end
            end

            if (trunc_flag && file_exist)  
                delete(filename); % check if delete succeeded?
                file_exist = false;
            end

            if (this.out_flag && ~file_exist) % create the file and add empty header
                this.fid = fopen(filename, 'w');
                if (this.fid<0), error(['[FileIO.FileIO]: Failed to create file "' filename '"...\n']); end
                this.writeHeader(); % write empty header
                fclose(this.fid);
            end

            % op_mode = std::ios::binary;
            if (this.in_flag), op_mode = 'r'; end
            if (this.out_flag), op_mode = 'r+'; end % add "r+" to avoid overriding previous contents

            this.fid = fopen(filename, op_mode);
            if (~this.fid), error(['[FileIO.FileIO]: Failed to open file "' filename '"... The file may not exist or it cannot be accessed...\n']); end

            this.readHeader();

        end
        
        %% Destructor.
        function delete(this)

            this.close();
            
        end
        
        %% Closes the file.
        %  If not called explicitly, it will be called eventually when this
        %  object goes out of scope.
        function close(this)
           
            if (this.fid > 0)
                fclose(this.fid);
                this.fid = -5;
            end
            
        end
        
        %% Prints the contents of the file on the console.
        function printHeader(this)
        
            name_len = 0;
            for k=1:length(this.name)
                if (length(this.name{k}) > name_len), name_len = length(this.name{k}); end
            end
            name_len = name_len + 5;

            type_name = cell(k,1);
            type_len = 0;
            for k=1:length(this.type) 
                type_name{k} = this.getFullTypeName(this.type{k}, this.sc_type{k});
                if (this.type{k}==FileIO.ARMA | this.type{k}==FileIO.EIGEN), type_name{k} = [type_name{k} this.getMatDim2str(k)]; end
                n = length(type_name{k});
                if (n > type_len), type_len = n; end
            end
            
            type_len = type_len + 5;

            horiz_line = repmat('-', 1, name_len + type_len + 6);

            fprintf('%s\n', horiz_line);
            fprintf('%s\n', ['Name' repmat(' ',1,name_len-length('Name')) 'Type' repmat(' ',1,type_len-length('Type')) 'Pos']);
            fprintf('%s\n', horiz_line);
            for k=1:length(this.name)
                fprintf('%s\n', [this.name{k} repmat(' ',1,name_len-length(this.name{k})) ...
                    type_name{k} repmat(' ',1,type_len-length(type_name{k})) num2str(this.i_pos{k})]);
            end
            
        end
            
        %% Read a variable from the file. 
        %  Openmode 'in' must be set.
        %  If the requested variable name is not found, an exception is thrown.
        %  @param[in] name_: the name of the variable.
        %  @param[out] s: the value of the variable.
        function s = read(this, name_)
            
            i = this.findNameIndex(name_);
            if (i<0), error(['[FileIO.read]: ' this.getErrMsg(FileIO.ENTRY_NOT_EXIST) ' "' name_ '"']); end

            t = this.type{i};
            matlab_type = this.getMatlabType(this.sc_type{i});
            fseek(this.fid, this.i_pos{i}, 'bof');

            if (t == FileIO.SCALAR)
                
                s = this.readScalar_(matlab_type);
                
            elseif (t == FileIO.ARMA | t == FileIO.EIGEN)
                
                [n_rows, n_cols] = this.readMatDim(name_);
                s = zeros(n_rows, n_cols, matlab_type);
                
                n_elem = n_rows*n_cols;
                buff = fread(this.fid, [n_elem, 1], [matlab_type '=>' matlab_type]);
                
                k = 1;
                for j=1:n_cols
                  for i=1:n_rows
                      s(i,j) = buff(k);
                      k = k + 1;
                  end
                end

            elseif (t == FileIO.CELL)
                
                [n_rows, ~] = this.readMatDim(name_);
                s = fread(this.fid, [n_rows 1], [matlab_type '=>' matlab_type]);
                % s = num2cell(s);
                
            elseif (t == FileIO.STRING)

                n_elem = fread(this.fid, 1, [this.long_t '=>' this.long_t]);
                s = fread(this.fid, [1 n_elem], 'char=>char');

            else
                
                error(['[FileIO.read]: Unsopported type "' matlab_type '".']);
                
            end
            
        end
        
        %% Write a variable to the file. 
        %  Openmode 'out' must be set.
        %  If the requested variable already exists in the file, an exception is thrown.
        %  @param[in] name_: the name of the variable.
        %  @param[in] s: the value of the variable.
        function write(this, name_, s)

            if (~this.out_flag), error(['[FileIO::write]: ' FileIO.getErrMsg(FileIO.INVALID_OP_FOR_OPENMODE) ': "' this.getOpenModeName() '"']); end
                
            if (this.findNameIndex(name_) >=0), error(['[FileIO::write]: ' this.getErrMsg(FileIO.DUPLICATE_ENTRY) ': "' name_ '"']); end
            
            if (isscalar(s)), this.writeScalar(name_, s);
            elseif (iscell(s)), this.writeCell(name_, s);
            elseif (ischar(s)), this.writeString(name_, s);
            elseif (sum(size(s)) > 2), this.writeMat(name_, s);
            else, error(['[FileIO.write]: Unsopported type "' class(s) '".']);
            end
            
        end
        
        %% Reads and returns all contents of the file as a struct.
        %  Openmode 'in' must be set.
        %  @param[out] s: a struct with attributes the names of the variables in the file.
        function s = readAll(this)
            
            s = struct();
            for i=1:length(this.name) s.(this.name{i}) = this.read(this.name{i}); end
            
        end
        
        
        %% Specify the c++ format for writing matrices.
        %  @param[in] mat_type: Type of c++ matrix \in {FileIO.ARMA, FileIO.EIGEN}
        function setMatType(this, mat_type)
            
            if (mat_type == FileIO.ARMA), this.MAT_TYPE = FileIO.ARMA;
            elseif (mat_type == FileIO.EIGEN), this.MAT_TYPE = FileIO.EIGEN;
            else error('[FileIO.setMatType]: Unsupported matrix type');
            end
            
        end
        
    end
 
    methods (Access = protected)
        
        function writeScalar(this, name_, s)

            k = length(this.name) + 1;
            
            this.name{k} = name_;
            this.type{k} = FileIO.SCALAR;
            this.sc_type{k} = this.findScalarType(class(s));
            this.name_map(name_) = k;
            fseek(this.fid, this.header_start, 'bof');
            this.i_pos{k} = ftell(this.fid);
            
            this.writeScalar_(s);
            
            this.header_start = ftell(this.fid);
            this.writeHeader(); % overwrites previous header
            
        end
        
        function writeCell(this, name_, v)

            [n_rows, n_cols] = size(v);
            
            if (n_rows~=1 & n_cols~=1)
                error('[FileIO::writeCell]: cell must be a vector.');
            end
            
            n_elem = max([n_rows, n_cols]);
            n_rows = cast(n_elem, this.long_t);
            n_cols = cast(1, this.long_t);

            t = this.findScalarType(class(v{1}));
            for i=1:n_elem
                t_i = this.findScalarType(class(v{i}));
                if (t ~= t_i)
                   error('[FileIO::writeCell]: cell entries must be of the same class(data type).');
                end
            end
            
            this.writeMatDim(class(v{1}), name_, FileIO.CELL, n_rows, n_cols);
            fwrite(this.fid, cell2mat(v), class(v{1}));

            this.header_start = ftell(this.fid);
            this.writeHeader(); % overwrites previous header
    
        end
        
        function writeString(this, name_, s)

            k = length(this.name) + 1;
            this.name{k} = name_;
            this.type{k} = FileIO.STRING;
            this.sc_type{k} = FileIO.ScType_NA;
            this.name_map(name_) = k;
            fseek(this.fid, this.header_start, 'bof');
            this.i_pos{k} = ftell(this.fid);
            
            n_elem = cast( length(s), this.long_t);
            
            fwrite(this.fid, n_elem, class(n_elem));
            fwrite(this.fid, s, class(s));

            this.header_start = ftell(this.fid);
            this.writeHeader(); % overwrites previous header
            
        end
        
        function writeMat(this, name_, m)
            
            [n_rows, n_cols] = size(m);
            n_rows = cast(n_rows, this.long_t);
            n_cols = cast(n_cols, this.long_t);
            
            this.writeMatDim(class(m), name_, FileIO.ARMA, n_rows, n_cols);
            
            buff = zeros(1, n_rows*n_cols, class(m));
            k=1;
            for j=1:n_cols
              for i=1:n_rows
                  buff(k) = m(i,j);
                  k = k + 1;
              end
            end
            fwrite(this.fid, buff, class(buff));

            this.header_start = ftell(this.fid);
            this.writeHeader(); % overwrites previous header
    
        end
        
        function writeMatDim(this, sc_t, name_, t, n_rows, n_cols)

            k = length(this.name) + 1;
            this.name{k} = name_;
            this.type{k} = t;
            this.sc_type{k} = this.findScalarType(sc_t);
            this.name_map(name_) = k;
            fseek(this.fid, this.header_start, 'bof');
            this.i_pos{k} = ftell(this.fid);

            this.writeScalar_(n_rows);
            this.writeScalar_(n_cols);
            
        end
        
        
        function writeHeader(this)
            
            fseek(this.fid, this.header_start, 'bof');

            i1 = ftell(this.fid);

            for k=1:length(this.name)
                name_len = length(this.name{k});
                fwrite(this.fid, name_len, this.long_t);
                fwrite(this.fid, this.name{k}, 'char');
                fwrite(this.fid, this.type{k}-this.type_offset, this.enum_t);
                fwrite(this.fid, this.sc_type{k}-this.scType_offset, this.enum_t);
                fwrite(this.fid, this.i_pos{k}, this.size_t);
            end

            i2 = ftell(this.fid);

            header_len = cast(i2-i1 + this.sizeof(this.long_t), this.long_t);
            fwrite(this.fid, header_len, class(header_len));
        
        end
        
        function readHeader(this)
            
            this.name_map = [];
            this.name = [];
            this.type = [];
            this.sc_type = [];
            this.i_pos = [];

            fseek(this.fid, 0, 'bof');
            i_start = ftell(this.fid);     
            fseek(this.fid, 0, 'eof');
            i_end = ftell(this.fid);
            
            if (i_start == i_end), error(['[FileIO.readHeader]: ' FileIO.getErrMsg(FileIO.EMPTY_HEADER)]); end

            header_len = cast(0, this.long_t);
            fseek(this.fid, -this.sizeof(header_len), 'eof');
            i_end = ftell(this.fid);
            header_len = fread(this.fid, 1, [class(header_len) '=>' class(header_len)]);
            if (header_len < 0), error(['FileIO.readHeader]: ' FileIO.getErrMsg(FileIO.CORRUPTED_HEADER)]); end

            this.name_map = containers.Map();
            
            fseek(this.fid, -header_len, 'eof');
            this.header_start = ftell(this.fid);
            while (ftell(this.fid) ~= i_end)
                len = fread(this.fid, 1, [this.long_t '=>' this.long_t]);
                name_i = fread(this.fid, [1 len], ['char=>char']);
                t = fread(this.fid, 1, [this.enum_t '=>' this.enum_t]) + this.type_offset;
                sc_t = fread(this.fid, 1, [this.enum_t '=>' this.enum_t]) + this.scType_offset;
                i = fread(this.fid, 1, [this.size_t '=>' this.size_t]);
                
                k = length(this.name)+1;
                this.type{k} = t;
                this.sc_type{k} = sc_t;
                this.i_pos{k} = i;
                this.name{k} = name_i;
                this.name_map(name_i) = k;
            end
            
        end

        function [n_rows, n_cols] = readMatDim(this, name_)

            i = this.findNameIndex(name_);
            if (i<0), error(['[FileIO.read]: ' this.getErrMsg(FileIO.ENTRY_NOT_EXIST) '"' name_ '"']); end

            fseek(this.fid, this.i_pos{i}, 'bof');

            n_rows = this.readScalar_(this.long_t);
            n_cols = this.readScalar_(this.long_t);

            if (n_rows<0 || n_cols<0)
              error(['[FileIO.read]: ' this.getErrMsg(FileIO.CORRUPTED_HEADER) ...
                  ': Negative matrix dimensions: "' name_ '" ' num2str(n_rows) ' x ' num2str(n_cols)] );
            end

        end
  
  
        function out = getMatDim2str(this, k)

            fseek(this.fid, this.i_pos{k}, 'bof');
            n_rows = this.readScalar_(this.long_t);
            n_cols = this.readScalar_(this.long_t);
            out = sprintf('(%i,%i)', n_rows, n_cols);

        end
           
        function s = readScalar_(this, type)

            s = fread(this.fid, 1, [type '=>' type]);

        end
        
        function writeScalar_(this, s)
  
            fwrite(this.fid, s, class(s));
            
        end
 
        function i = findNameIndex(this, name_)
        
            if (this.name_map.isKey(name_)), i = this.name_map(name_);
            else, i = -1;
            end
                
        end
        
        function t = findScalarType(this, sc_t)
  
            if (strcmpi(sc_t,'int8')), t = FileIO.BOOL;
            elseif (strcmpi(sc_t,'int32')), t = FileIO.INT;
            elseif (strcmpi(sc_t,'uint32')), t = FileIO.UINT;
            elseif (strcmpi(sc_t,'int32')), t = FileIO.LONG;
            elseif (strcmpi(sc_t,'uint32')), t = FileIO.ULONG;
            elseif (strcmpi(sc_t,'int64')), t = FileIO.LLONG;
            elseif (strcmpi(sc_t,'uint64')), t = FileIO.ULLONG;
            elseif (strcmpi(sc_t,'single')), t = FileIO.FLOAT;
            elseif (strcmpi(sc_t,'double')), t = FileIO.DOUBLE;
            else, t = FileIO.UNKNOWN;
            end
    
            if (t == FileIO.UNKNOWN), error(['[FileIO.findScalarType]: ' this.getErrMsg(FileIO.UNKNOWN_TYPE)]); end
            
        end
          
        function op_mode = getOpenModeName(this)
            
            if (this.in_flag), op_mode = 'in';
            elseif (this.out_flag), op_mode = 'out';
            else, op_mode = 'in|out';
            end
  
        end
    end
    
    methods (Static, Access = protected)
        
        function err_msg = getErrMsg(err_id)
            
            if (err_id == FileIO.EMPTY_HEADER)
                err_msg = 'EMPTY_HEADER';
            elseif (err_id == FileIO.CORRUPTED_HEADER)
                err_msg = 'CORRUPTED_HEADER';
            elseif (err_id == FileIO.ENTRY_NOT_EXIST)
                err_msg = 'ENTRY_NOT_EXIST';
            elseif (err_id == FileIO.DUPLICATE_ENTRY)
                err_msg = 'DUPLICATE_ENTRY';
            elseif (err_id == FileIO.TYPE_MISMATCH)
                err_msg = 'TYPE_MISMATCH';
            elseif (err_id == FileIO.DIMENSIONS_MISMATCH)
                err_msg = 'DIMENSIONS_MISMATCH';
            elseif (err_id == FileIO.UNKNOWN_TYPE)
                err_msg = 'UNKNOWN_TYPE';
            end

        end
        
        
        function sc_t_name = getMatlabType(sc_t)
            
            if (sc_t == FileIO.BOOL)
                sc_t_name = 'int8';
            elseif (sc_t == FileIO.INT)
                sc_t_name = 'int32';
            elseif (sc_t == FileIO.UINT)
                sc_t_name = 'uint32';
            elseif (sc_t == FileIO.LONG)
                sc_t_name = 'int32';
            elseif (sc_t == FileIO.ULONG)
                sc_t_name = 'uint32';
            elseif (sc_t == FileIO.LLONG)
                sc_t_name = 'int64';
            elseif (sc_t == FileIO.ULLONG)
                sc_t_name = 'uint64';
            elseif (sc_t == FileIO.FLOAT)
                sc_t_name = 'single';
            elseif (sc_t == FileIO.DOUBLE)
                sc_t_name = 'double';
            elseif (sc_t == FileIO.ScType_NA)
                sc_t_name = '';
            else
                sc_t_name = 'N/A';
            end

        end
          
        
        function sc_t_name = getScalarTypeName(sc_t)
            
            if (sc_t == FileIO.BOOL)
                sc_t_name = 'int8';
            elseif (sc_t == FileIO.INT)
                sc_t_name = 'int';
            elseif (sc_t == FileIO.UINT)
                sc_t_name = 'uint';
            elseif (sc_t == FileIO.LONG)
                sc_t_name = 'long';
            elseif (sc_t == FileIO.ULONG)
                sc_t_name = 'ulong';
            elseif (sc_t == FileIO.LLONG)
                sc_t_name = 'int64';
            elseif (sc_t == FileIO.ULLONG)
                sc_t_name = 'uint64';
            elseif (sc_t == FileIO.FLOAT)
                sc_t_name = 'float';
            elseif (sc_t == FileIO.DOUBLE)
                sc_t_name = 'double';
            elseif (sc_t == FileIO.ScType_NA)
                sc_t_name = '';
            else
                sc_t_name = 'N/A';
            end

        end

        function t_name = getTypeName(t)
            
            if (t == FileIO.SCALAR)
                t_name = 'scalar';
            elseif (t == FileIO.ARMA)
                t_name = 'mat';
            elseif (t == FileIO.EIGEN)
                t_name = 'mat';
            elseif (t == FileIO.CELL)
                t_name = 'cell';
            elseif (t == FileIO.STRING)
                t_name = 'string';
            else
                t_name = 'N/A';
            end

        end
        
        function full_t_name = getFullTypeName(type, sc_type)
            
            sc_t = FileIO.getScalarTypeName(sc_type);
            if (isempty(sc_t)), full_t_name = FileIO.getTypeName(type);
            else, full_t_name = [FileIO.getTypeName(type) '<' sc_t '>'];
            end
            
        end
        
        function nbytes = sizeof(var)

            narginchk(1, 1);

            if (ischar(var))
                precision = var;
                try
                    z = zeros(1, precision);
                catch
                    error('Unsupported class for finding size');
                end
            else
               z = var;
            end

            w = whos('z');
            nbytes = w.bytes;

        end

    end

    properties (Access = protected)

        fid % the stream for reading/writing to the file

        header_start % the position where the header starts in the file

        name_map % maps the names of the data to their index in the following vectors
        name % cell(vector) with the names of the data contained in the file
        type % cell(vector) with the type (see @Type) of the data contained in the file
        sc_type % cell(vector) with the scalar type (see @ScalarType) of the data contained in the file
        i_pos % cell(vector) with the position of the data in the file
        
        MAT_TYPE % the matrix type (AMRA, EIGEN) assigned to matlab matrices when 'write' is called

        % static const char* error_msg[];
        % static const char *TypeName[];
        % static const char *ScalarTypeName[];
        
        in_flag; % true when the open_mode is "in"
        out_flag; % true when the open_mode is "out"

    end
    
    properties (Constant, Access = public)
        
        % enum OpenMode
        in = bitshift(1, 0);     % Allow input operations.
        out = bitshift(1, 1);    % Allow output operations.
        trunc = bitshift(1, 2);   % Discard all previous contents from the file.
        % binary = 1<<3

        % enum Type
        SCALAR = int32(1001);    % scalar, see @ScalarType
        ARMA = int32(1002);      % arma (Mat, Col, Row), Eigen (Matrix, Vecotr, RowVector)
        EIGEN = int32(1003);
        CELL = int32(1004);      % std::vector
        STRING = int32(1005);      % std::string
        % CUSTOM ?


        % enum ScalarType
        BOOL = int32(2001);      % int8 (bool)
        INT = int32(2002);       % int32 (int)
        UINT = int32(2003);      % uint32 (unsigned int)
        LONG = int32(2004);      % int32 (long)
        ULONG = int32(2005);     % uint32 (unsigned long)
        LLONG = int32(2006);     % int64 (long long)
        ULLONG = int32(2007);    % uint64 (unsigned long long)
        FLOAT = int32(2008);     % single (float)
        DOUBLE = int32(2009);    % double (double)
        UNKNOWN = int32(2010);   % unsupported type
        ScType_NA = int32(2011);   % for Types where the scalar type is redundant, like strings

      
        % enum Error
        EMPTY_HEADER = int32(100);
        CORRUPTED_HEADER = int32(101);
        ENTRY_NOT_EXIST = int32(102);
        DUPLICATE_ENTRY = int32(103);
        TYPE_MISMATCH = int32(104);
        DIMENSIONS_MISMATCH = int32(105);
        UNKNOWN_TYPE = int32(106);
        
        long_t = 'int64';
        size_t = 'uint64';
        enum_t = 'int32'
        
        type_offset = int32(1001);
        scType_offset = int32(2001);
        
    end
    
end