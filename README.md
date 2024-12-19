# Coherent Lagrangian Vortices

## Background

A script for computing *coherent Lagrangian vortex* (CLV) boundaries from flows defined on curved surfaces in arbitrary coordinates.  It uses geodesic vortex detection[^1] as implemented in the Julia package [CoherentStructures.jl](https://github.com/CoherentStructures/CoherentStructures.jl), a software tool for automated CLV detection based on line-field theory[^2] [^3].  The script `clv.jl` allows one to define the metric consistent with coordinate representation of the velocity field data input.  It also fixes some issues with interpolations as originally implemented in the original [CoherentStructures.jl](https://github.com/CoherentStructures/CoherentStructures.jl).  The input and output formats can be modified.  While the data included for testing pertain to stratospheric flow, the surfaces on which flow is defined can include coastlines.

## Usage

### 1.

Install [Julia](https://julialang.org/). For most users, it is enough to run the following command in a terminal. See [here](https://github.com/JuliaLang/juliaup) for installation instructions if this is not sufficient.

On Mac/Linux:

```sh
curl -fsSL https://install.julialang.org | sh
```

On Windows:

```sh
winget install julia -s msstore
```

### 2. 

Download this repository as either a [zipped folder](https://github.com/70Gage70/CoherentLagrangianVortices/archive/refs/heads/master.zip), or by [cloning](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) it and `cd` to that folder.

### 3.

From the shell, run the following to install all the required dependencies:

```sh
julia init.jl
```

From the shell, run the following to execute the main calculations on the provided test data in `ERA5-600K-May.nc`:

```sh
julia --project=. --threads=auto clv.jl
```

This should take roughly 1 minute and generate `output.png` showing the location of singularities and the outline of the barrier as well as `output.nc` which has the actual data.

### Customization

The data can be output to `.mat` by commenting out the indicated sections of `clv.jl`. To use custom data, modify the `fileIn` variable to the name of your file.

## License

This project is licensed under the GNU General Public License v3.0. See the [LICENSE](LICENSE.md) file for details.


# References

[^1]: Haller, G., & Beron-Vera, F. J. (2013). Coherent Lagrangian vortices: The black holes of turbulence. Journal of fluid mechanics, 731, R4.

[^2]: Karrasch, D., Huhn, F., & Haller, G. (2015). Automated detection of coherent Lagrangian vortices in two-dimensional unsteady flows. Proceedings of the Royal Society A: Mathematical, Physical and Engineering Sciences, 471(2173), 20140639.

[^3]: Karrasch, D., & Schilling, N. (2020). Fast and robust computation of coherent Lagrangian vortices on very large two-dimensional domains. The SMAI journal of computational mathematics, 6, 101-124.

