# Jlab Simulation Package for Electron Cooling

## About JSPEC
JSPEC is an open source C++ package for numerical simulations on the electron cooling process, including the intrabeam scattering (IBS) effect, developed at [Jefferson Lab (JLab)](http://www.jlab.org). 



THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## How to compile
#### Compile using code blocks IDE

JSPEC is developed using the code blocks IDE. If you are using the same IDE, just open the project file "jspec.cbp" in the cbp folder, select a C++ compiler and build it. Make sure your compiler supports C++11 standard and the option std=c++11 is on. 

#### Compile using the makefile

A makefile is provided (tested in Ubuntu 18.04). In the project folder, run the following commands:

> `make`

The executable file called "jspec" will be generated.

## How to run

To run JSPEC, you can put your input file in the same folder with the JSPEC executable file and run the following commands in the folder:

> ` jspec.exe inputfilename` 

You also need another file in [*MAD X*](https://madx.web.cern.ch/madx/) tfs format, which defines the ion ring optics. You can put it in the same folder too. About how to write your input file, please refer to the JSPEC User Manual. 

Besides running the program in the command line as aforementioned, a windows batch file "jspec-dragfile.bat" is also provided. Putting this bat file together with the executable jspec file "jspec.exe", the math parser dynamic library file "muparser.dll", the tfs format lattice file and the input script file in the same folder, one can drag the  input script file onto the batch file to run the computation/simulation. 

## Acknowledgement

Authors of [**BETACOOL**](http://betacool.jinr.ru/), we learned a lot from BETACOOL. 

Authors of [**muParser**](http://beltoforion.de/article.php?a=muparser),  which is used in building the text-based UI. 

Dr. David Bruwihler, Dr. Paul Moeller and Dr. Stephen Coleman at [*Radiasoft*](http://radiasoft.net/), who developed an [*online version of JSPEC with GUI*](https://beta.sirepo.com/#/jspec) on their cloud server, [*Sirepo*](https://beta.sirepo.com/). 

My colleagues at Jefferson Lab. 



## Contact the authors 

Dr. He Zhang at [*Jefferson Lab*](www.jlab.org) by hezhang.AT.jlab.org. 