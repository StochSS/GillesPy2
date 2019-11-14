Installation
############

GillesPy2 can be installed on your computer using different methods.


Preferred method: using PyPI
****************************

Using Python 3 on **Linux**, **macOS**, and **Windows** operating systems, you should be able to install GillesPy2 using the package management system `pip <https://pip.pypa.io/en/stable/installing>`_ by typing the following commands in a command shell interpreter::

    python3 -m pip install gillespy2 --user --upgrade


Alternative methods: using the code repository
**********************************************

As an alternative to getting it from PyPI, you can instruct ``pip`` to install GillesPy2 directly from the `GitHub repository for GillesPy2 <https://github.com/gillesPy2/GillesPy2>`_ by using the following command::

    python3 -m pip install git+https@github.com:GillesPy2/GillesPy2.git --user --upgrade


As a final alternative, you can first use ``git`` to clone a copy of the GillesPy2 source tree from the GitHub repository and then install it using that copy::

    git clone --recursive https@github.com:GillesPy2/GillesPy2.git
    cd GillesPy2
    python3 -m pip install . --user --upgrade
