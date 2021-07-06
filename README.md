# Rustpond
Rustpond -- A Rust port of Nanopond 1.9, a teeny tiny artificial life virtual machine\
Rust code by Josh Franklin\
Original C implementation of Nanopond by Adam Ierymenko (https://github.com/adamierymenko/nanopond) \
(v1.9 link: https://web.archive.org/web/20170816191509/http://adam.ierymenko.name/nanopond.shtml)

Rustpond is an implementation of Nanopond 1.9. 

## Installation
Requires a Rust 2018 compatible installation of cargo and rustc. Instructions for downloading and installing Rust are available at https://www.rust-lang.org/learn/get-started. You will need to follow the link for "Other installation methods" if you are using Windows.

```bash
git clone https://github.com/joshlfrankl/rustpond.git
cargo build --release
```
The rustpond executable will be in the /target/release directory. If you do not wish to use git, download this repository as a zip file, then enter the cargo command inside the unzipped directory.

## Usage
```bash
./rustpond --help

Rustpond 0.9.0
Josh Franklin <joshlfrankl@gmail.com>
A Rust port of Nanopond, a tiny artificial life VM. Original Nanopond C implementation written by Adam Ierymenko.

USAGE:
    rustpond [FLAGS] [OPTIONS]

FLAGS:
    -b, --benchmark     Run a benchmark with a consistent seed for one billion updates; ignores all other flags
    -h, --help          Prints help information
    -V, --version       Prints version information
    -p, --write_pngs    Writes PNG files every REPORT_FREQUENCY updates; default is false

OPTIONS:
    -o, --out_path <out_path>    Output filepath to save PNGs; default is the current directory
    -s, --seed <seed>            Random seed to use; default is a randomly generated seed
    -u, --updates <updates>      Number of updates to run; default is 2^64 - 1
```
Use the benchmark flag for performance testing. The results should be comparable across all OS and hardware configurations.

If the -p flag is passed to Rustpond, it will write images of the pond to sequentially numbered png files. These files can be converted into a video using [ffmpeg](http://ffmpeg.org/), e.g.:
```bash
ffmpeg -r 20 -f image2 -s 512x512 -i rp_%d.png -vcodec libx264 -crf 16 rp_movie.mp4
```
You can control the location these files are written to using the -o flag.

The -s flag is used to specify a seed for reproducibility. The seed should be an integer between 0 and 2^64 - 1. For a given seed, Rustpond should emit identical output each time it is run.

Rustpond will always log information about the simulation to a csv file titled "output.csv". You can change the location of the output file using the -o flag.

## Configuration
As in the original Nanopond implementation, most of the simulation parameters are constants at the top of the main.rs file. This allows for the compiler to make some optimizations at the cost of reduced flexibility. In the future, I may alter Rustpond to allow command line or configuration file control over the simulation parameters. Currently, you will have to re-build Rustpond (cargo build --release) after any changes to the parameters. Be sure to use the updated executable after re-building.

## License
[GNU GPLv3](http://www.gnu.org/licenses/gpl-3.0.txt)






