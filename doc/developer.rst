Developing on GBasis
====================

Modifying GBasis is somewhat possible. It requires Python, Cython, and C++
familiarity. Because of the mixed-language nature of the code, development
tools don't work comprehensively. You can still use Python dev tools for
Python code, but it is unable to analyze the sections which jump to C++.

Since it is not possible to step through ``GBasis`` with a debugger, we
made diagrams explaining the rough code execution when calling
:code:`GOBasis.compute_overlap` and
:code:`GOBasis.compute_electron_repulsion`.

First, we cover how to call C++ from Python.

.. _calling_cxx:

Calling C++ from Python using Cython
------------------------------------
*Skip this if you are intimately familiar with Cython*

We give a quick primer to Cython in this section (this is in no way a
replacement for the Cython documentation!).

Developers use Cython to accelerate Python by writing ``.pyx`` files.
These files are compiled by Cython into C++ files, which generate a shared
object that CPython can import. If you wish to run your own C++ code in
Python, you must import it into the ``.pyx`` file.

Importing C++ code into a ``.pyx`` is possible with the help of a ``.pxd``
file. Think of it as the Cython equivalent of the C++ ``.h`` header.
Once a C++ object has been declared in a ``.pxd`` header, you can call
it from within the ``.pyx``.

So for example if we wanted to expose the :code:`compute_grid_point1`
function from the C++ class :code:`GBasis` in the header ``gbasis.h``:

.. code-block:: c++

    class GBasis {
        void compute_grid_point1(double *output, double *point, GB1GridFn *grid_fn);
    };

We would generate a ``gbasis.pxd``

.. code-block:: python

    cdef extern from "gbasis.h":
        cdef cppclass GBasis:
            void compute_grid_point1(double* output, double* point, fns.GB1DMGridFn* grid_fn)

Subsequently, we could call it from the ``.pyx`` file

.. code-block:: python

    cimport gbasis

    cdef class GBasis:
        cdef gbasis.GBasis* _this

        def grid_point1(self, double[::1] output not None,
                                double[::1] point not None,
                                GB1DMGridFn grid_fn not None):
            # ... type safety checks go here
            self._this.compute_grid_point1(&output[0], &point[0], grid_fn._this)

We could then call it from Python

.. code-block:: python

    obasis = GBasis()
    obasis.grid_point1(...)

Note that the Cython syntax is a mix of C++ and Python. You are responsible for handling some
pointer dereferencing yourself. Please refer to the
`Cython docs <https://cython.readthedocs.io/en/latest/>`_ for more details.

.. warning::

    The type safety checks within ``.pyx`` files are *important*. If you pass incompatible
    objects to C++ code, you will get a segfault at runtime, which is extremely hard to debug
    since gdb and pdb won't work. You won't even get a traceback. The compiler will also not
    warn you when compiling the Cython code. Pay extra close attention when writing that code!

GBasis structure
----------------

The previous section gives the general idea of how code flows within GBasis. All the user
interfaces are from within Python, which interprets compiled Cython code, which calls C++
classes/functions. On a user's machine, the C++ and Cython code is compiled into a shared
object file and CPython runs it as a binary. Another common use case in GBasis is to expose
C++ code for unit testing with ``nosetests``. This practice is a relic from the original
HORTON code and should not be emulated if possible. We now use google-test for direct
C++ unit testing. Regardless, we document the code so developers can modify it.

This is an overview of how code flows when calculating integrals:

.. code-block:: python

    from gbasis import get_gobasis

    # ... get coordinates, numbers, basis
    obasis = get_gobasis(coordinates, number, basis)
    olp = obasis.compute_overlap()
    er = obasis.compute_electron_repulsion()

First the graph for generating the GOBasis instance:

.. code-block:: python

    obasis = get_gobasis(coordinates, number, basis)

.. mermaid::

    graph TD

    subgraph Python
    c(coordinates) --> A
    n("atomic numbers") --> A
    b("basis set") -->A
    A["get_gobasis (gobasis.py)"]
    end
    subgraph Cython

    A --> C("GOBasis (cext.pyx)")
    B("GBasis (cext.pyx)") -.- C
    end

    subgraph C++
    gbasis("GBasis (gbasis.cpp)") -.- D
    C --> D("GOBasis (gbasis.cpp)")
    end

    subgraph Legend
    cl("Class (filename.py)") -->|code flow| fn[function]
    par(Parent) -.-|inheritance| ch(Child)
    end

Implementing 1-electron integrals
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After generating the GOBasis instance, you can ask for an overlap integral:

.. code-block:: python

    olp = obasis.compute_overlap()

.. mermaid::

    graph TD

    START[" "] --> olp
    subgraph GOBasis : Gbasis  gbasis.cpp
        olp[compute_overlap]
    end

    subgraph GBasis gbasis.cpp
        cti[compute_two_index]
    end

    subgraph GB2Integral ints.cpp
        ctp[cart_to_pure]
        reset[reset]
    end

    subgraph  IterGB2 iter_gb.cpp
        upsh[update_shell]
        store
        incsh[inc_shell]
    end

    loop{"loop over primitives"}

    olp --> cti
    cti --> upsh
    upsh --> reset
    reset --> loop
    loop --> ctp
    ctp --> store
    store -->|shells remain| incsh
    store -->|no shells left| END[" "]
    incsh --> reset

    style START fill:#FFFFFF, stroke:#FFFFFF;
    style END fill:#FFFFFF, stroke:#FFFFFF;

The loop over the primitives is as follows:

.. mermaid ::

    graph LR

    subgraph GB2OverlapIntegral:GB2Integral ints.cpp
        add[add]
    end

    subgraph  IterGB2 iter_gb.cpp
        uppr[update_prim]
        incpr[inc_prim]
    end

    START[" "] --> uppr
    uppr --> add
    add --> incpr
    incpr -->|primitives remain| add
    add -->|no primitives left|END[" "]

    style START fill:#FFFFFF, stroke:#FFFFFF;
    style END fill:#FFFFFF, stroke:#FFFFFF;

Astute readers will note that the only section which references the overlap integral explicitly
is in the :code:`add` function within GB2OverlapIntegral in the primitives loop. Other
1-electron integrals have similar structure. Thanks to the object oriented nature of this
library, it is very easy to replace it.

The following is the entirety of the GB2OverlapIntegral class header:

.. code-block:: python

    class GB2OverlapIntegral : public GB2Integral {
     public:
      explicit GB2OverlapIntegral(long max_shell_type) : GB2Integral(max_shell_type) {}

      virtual void add(double coeff, double alpha0, double alpha1, const double *scales0, const double *scales1);
    };

There is only one function implemented, the overlap integral kernel within :code:`add`,
and the rest is inherited from the parent GB2Integral class. **To implement your own 1-electron
integrals, all you need to do is to make a new child class of GB2Integral and implement the
add function.** You must also expose said function in ``ints.pxd`` and ``cext.pyx`` as mentioned
in the section :ref:`calling_cxx`. A complete implementation which
will be acceptable for merging into master will also have unit tests in ``test_ints.py`` and follow
coding standards. More details on merging contributions are available in :ref:`dev_building`

Implementing 2-electron integrals
---------------------------------

