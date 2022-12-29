io_
=========

.. image:: docs/_static/example.png
   :alt: io_ in action
   :align: left

.. -*- inclusion-marker-for-sphinx-docs -*-

io__ is a header-only C++ library for printing colored messages to the
terminal. Written just for fun with a help of `the Force`_. io_ uses
`ANSI color formatting`_, so you can use it on every system that is used such
terminals (most \*nix systems, including Linux and Mac OS). On Windows, WinAPI
is used instead but some limitations are applied.

It's licensed under the BSD (3-clause) License. That basically means:
do whatever you want as long as copyright sticks around.

.. _io_: https://github.com/ikalnitsky/io_
.. _the Force: http://starwars.wikia.com/wiki/The_Force
.. _ANSI color formatting: http://en.wikipedia.org/wiki/ANSI_escape_code#Colors


Installation
------------

Add ``io_.hpp`` to the project and use provided stream manipulators
from the ``io_`` namespace.


How to use?
-----------

It's very easy to use. The idea is based on the use of C++ stream
manipulators. The typical «Hello World» application is below:

.. code:: c++

    #include <iostream>
    #include <io_/io_.hpp>

    int main(int /*argc*/, char** /*argv*/)
    {
        std::cout << io_::red << "Hello, Colorful World!" << std::endl;
        return 0;
    }

The application above prints a string with red. It's obvious, isn't it?
There is a one problem that is not obvious for the unexperienced users.
If you write something this way:

.. code:: c++

    std::cout << io_::red << "Hello, Colorful World!" << std::endl;
    std::cout << "Here I'm!" << std::endl;

the phrase «Here I'm» will be printed with red too. Why? Because you don't
reset io_'s setting. So if you want to print text wit default terminal
setting you have to reset io_'s settings. It can be done by using
``io_::reset`` manipulator:

.. code:: c++

    std::cout << io_::red << "Hello, Colorful World!" << std::endl;
    std::cout << io_::reset << "Here I'm!" << std::endl;

By default, io_ ignores any colors for non-tty streams
(e.g. ``std::stringstream``), so:

.. code:: c++

    std::stringstream ss;
    ss << io_::red << "unicorn";
    std::cout << ss.str();

would print «unicorn» using default color, not red. In order to change this
behaviour one can use ``io_::colorize`` manipulator that enforce colors
no matter what.


What manipulators are supported?
--------------------------------

The manipulators are divided into four groups:

* *foreground*, which changes text color;
* *background*, which changes text background color;
* *attributes*, which changes some text style (bold, underline, etc);
* *control*, which changes io_'s behaviour.


Foreground manipulators
.......................

#. ``io_::grey``
#. ``io_::red``
#. ``io_::green``
#. ``io_::yellow``
#. ``io_::blue``
#. ``io_::magenta``
#. ``io_::cyan``
#. ``io_::white``


Background manipulators
.......................

#. ``io_::on_grey``
#. ``io_::on_red``
#. ``io_::on_green``
#. ``io_::on_yellow``
#. ``io_::on_blue``
#. ``io_::on_magenta``
#. ``io_::on_cyan``
#. ``io_::on_white``


Attribute manipulators
......................

(so far they aren't supported on Windows)

#. ``io_::bold``
#. ``io_::dark``
#. ``io_::underline``
#. ``io_::blink``
#. ``io_::reverse``
#. ``io_::concealed``

Control manipulators
....................

(so far they aren't supported on Windows)

#. ``io_::colorize``
#. ``io_::nocolorize``
