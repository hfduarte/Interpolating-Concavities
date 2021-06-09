# Interpolating-Concavities

An algorithm to interpolate concavities.

J. Duarte and M. McKenney. Interpolating Concavities. ACM SIGAPP Applied Computing Review. vol 21, issue 3, 2021.

/gui contains a Java GUI for visulization.

![GUI](https://github.com/hfduarte/Interpolating-Concavities/blob/main/gui.png)

/src contains the source code for the algorithm to interpolate concavities as presented in the paper mentioned. 
The algorithm was implemented using python 3. Because of a dependency on the Scientific Python Geometric Algorithms Library we use miniconda3. 
lcip_opt.py contains the implementation described on the paper. 
lcip_nopt.py contains an older implementation.

/dep contains the GUI dependencies.

/datasets contains the datasets used for testing.

/exp contains images of some experiments.

## Installation

### GUI Installation

```bash
cd /home/user
git clone https://github.com/hfduarte/Interpolating-Concavities.git
```

Open eclipse choose File > Import.

In the Import Window select General > File System and click Next.

In the next window in From Directory input /home/user/Interpolating-Concavities/gui.

Choose Select All and click Finish.

### Source Installation

Install miniconda3

Install the dependencies, shapely and skgeom, using miniconda.

The tests can be runned using the GUI or directly by calling the python script as:

```bash
python lcip_opt.py <number_of_tests>
```

where number_of_tests is the number of times a test case is runned (tested).

## Dependencies

GUI Dependencies

[JTS Topology Suite](https://github.com/locationtech/jts)

The JTS Topology Suite is a Java library for creating and manipulating vector geometry. It also provides a comprehensive set of geometry test cases, and the TestBuilder GUI application for working with and visualizing geometry and JTS functions.

Source Code Dependencies

[Miniconda3](https://docs.conda.io/en/latest/miniconda.html)

[Shapely](https://pypi.org/project/Shapely/)

[skgeom](https://pythonawesome.com/scientific-python-geometric-algorithms-library/)

The Scientific Python Geometric Algorithms Library.

## Authors

José Duarte: hfduarte@ua.pt

Mark McKenney

## License

[MIT](https://choosealicense.com/licenses/mit/)
