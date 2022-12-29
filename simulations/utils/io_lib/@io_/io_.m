%% FileIO class
% 


classdef io_

    methods (Static, Access = public)
        
        %  \brief Reads a 2D matrix from stream 'fid'.
        %  \details Reads first the number of rows and columns and then the elements of the matrix row by row.
        %  @param[in] fid: The input stream.
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] type: The class type of 'm' (optional, default = 'double').
        %  @param[out] m: The 2D matrix
        function m = read_mat(fid, binary, type)

            if (nargin < 2), binary = false; end % text format
            if (nargin < 3), type = 'double'; end

            n_rows = io_.read_scalar(fid, binary, 'int64');
            n_cols = io_.read_scalar(fid, binary, 'int64');

            m = io_.read_mat_(fid, n_rows, n_cols, binary, type);

        end
        
        %  \brief Reads a scalar value from stream \a in.
        %  \details Reads a scalar value in the format specified by \a binary flag.
        %  @param[out] scalar scalar value.
        %  @param[in] fid: The input stream.
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] type: The class type of 'scalar' (optional, default = 'double').
        function scalar = read_scalar(fid, binary, type)

            if (nargin < 2), binary = false; end % text format
            if (nargin < 3), type = 'double'; end

            if (binary)
                scalar = fread(fid, 1, type);
            else     
                s = fscanf(fid,'%s', 1);
                scalar = str2num(s);
            end

        end

        %  \brief Reads a vector of 2D matrices from stream 'fid'
        %  \details Reads the number of 2D matrices, the number of rows and cols and the elements of each 2D matrix row by row.
        %  @param[in] fid: The input stream.
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] type: The class type of 'm' (optional, default = 'double').
        %  @param[out] m: cell array where each cell has a 2D matrix
        function m = read_vec_mat(fid, binary, type)

            if (nargin < 2), binary = false; end % text format
            if (nargin < 3), type = 'double'; end

            n_mat = io_.read_scalar(fid, binary, 'int64');

            m = cell(n_mat,1);

            for k=1:n_mat
                m{k} = io_.read_mat(fid, binary, type);
            end

        end
        
        %  \brief Writes a 2D matrix in stream 'fid'.
        %  \details Writes first the number of rows and columns and then the elements of the matrix row by row.
        %  @param[in] m: The 2D matrix
        %  @param[in] n_rows: The number of rows
        %  @param[in] n_cols: The number of columns
        %  @param[in] fid: The output stream (optional, default = 1 for output to screen).
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] precision: Precision in txt format (optional, default = 6).
        function write_mat(m, fid, binary, precision)

            if (nargin < 2), fid = 1; end % write to screen
            if (nargin < 3), binary = true; end % text format
            if (nargin < 4), precision = 6; end % used in text format

            n_rows = int64(size(m,1));
            n_cols = int64(size(m,2));

            io_.write_scalar(n_rows, fid, binary, precision);
            if (~binary), fprintf(fid, '\n'); end
            io_.write_scalar(n_cols, fid, binary, precision);
            if (~binary), fprintf(fid, '\n'); end

            io_.write_mat_(m, n_rows, n_cols, fid, binary, precision);

        end
        
        %  \brief Writes a scalar value in stream 'fid'.
        %  \details Writes the scalar value in the format specified by 'binary' flag. If the format is binary,
        %           the class type of 'scalar' is used to determine the number of bits to use.
        %  @param[in] scalar: scalar value
        %  @param[in] fid: The output stream (optional, default = 1 for output to screen).
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] precision precision in txt format (optional, default = 6).
        function write_scalar(scalar, fid, binary, precision)

            if (nargin < 2), fid = 1; end % write to screen
            if (nargin < 3), binary = false; end % text format
            if (nargin < 4), precision = 6; end % used in text format

            if (binary)
                fwrite(fid, scalar, class(scalar));
            else     
                s = num2str(scalar, precision);
                fprintf(fid, '%s', s);
            end

        end

        %  \brief Writes a vector of 2D matrices in stream 'fid'
        %  \details Writes the number of 2D matrices, and thern the number of rows and cols and the elements of each 2D matrix row by row.
        %  @param[out] m: cell array where each cell has a 2D matrix
        %  @param[in] fid: The output stream (optional, default = 1 for output to screen).
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] precision precision in txt format (optional, default = 6).
        function write_vec_mat(m, fid, binary, precision)

            if (nargin < 2), fid = 1; end % write to screen
            if (nargin < 3), binary = false; end % text format
            if (nargin < 4), precision = 6; end % used in text format

            n_mat = int64(length(m));

            io_.write_scalar(n_mat, fid, binary, precision);
            if (~binary), fprintf(fid, '\n'); end

            for k=1:n_mat
                io_.write_mat(m{k}, fid, binary, precision);
            end
        end

    end
    
    methods (Static, Access = protected)
        
        %  \brief Reads a 2D matrix from stream 'fid'.
        %  \details Reads the elements of the matrix row by row.
        %  @param[in] fid: The input stream.
        %  @param[in] n_rows: The number of rows.
        %  @param[in] n_cols: The number of columns.
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] type: The class type of 'm' (optional, default = 'double').
        %  @param[out] m: The 2D matrix
        function m = read_mat_(fid, n_rows, n_cols, binary, type)

            if (nargin < 4), binary = false; end % text format
            if (type < 5), type = 'double'; end

            m = zeros(int64(n_rows),int64(n_cols));

            if (binary)    
                for i=1:n_rows
                    m(i,:) = fread(fid, [1 n_cols], type);
                end  
            else
                for i=1:n_rows
                    m(i,:) = fscanf(fid,'%f', n_cols);
                end

            end

        end
        
        %  \brief Writes a 2D matrix in stream 'fid'.
        %  \details Writes first the number of rows and columns and then the elements of the matrix row by row.
        %  @param[in] m: The 2D matrix
        %  @param[in] n_rows: The number of rows
        %  @param[in] n_cols: The number of columns
        %  @param[in] fid: The output stream (optional, default = 1 for output to screen).
        %  @param[in] binary: Flag indicating the format (true for binary, false for text, optional, default = false).
        %  @param[in] precision: Precision in txt format (optional, default = 6).
        function write_mat_(m, n_rows, n_cols, fid, binary, precision)

            if (nargin < 4), fid = 1; end % write to screen
            if (nargin < 5), binary = false; end % text format
            if (nargin < 6), precision = 6; end % used in text format

            if (binary)
                for i=1:n_rows
                    fwrite(fid, m(i,:), class(m));
                end 
            else 
                for i=1:n_rows
                    s = num2str(m(i,:), precision);
                    fprintf(fid, '%s\n', s);
                end
            end

        end

        
    end
 
end
