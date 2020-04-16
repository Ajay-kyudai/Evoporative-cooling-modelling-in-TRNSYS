# Evoporative-cooling-modelling-in-TRNSYS
Modelling evaporative heat transfer by combining TRNSYS and MATLAB
* well known building energy simulation tool Trnsys can't model evaporative cooling from building surface due to moisture trasfer from building envelop.
* This project devises a new method to enble Trnsys to evaluate evaporative cooling. Wet roof with evaporation is modelled in Matlab and building is modelled in Trnsys. Then both are combined using a method shown in the figures.
## Pictorial Explanation of the alorithm
![](Algorithm%20explanation/Algorithm%20explanation.jpg)
![](Algorithm%20explanation/Algorithm%20explantion2.jpg) ![](Algorithm%20explanation/Algorithm%20explantion3.jpg)
![](Algorithm%20explanation/TRNSYS%20connection%20diagram.jpg)
## How to run
* clone or download the files and save it in your Trnsys repository
* As it contains also Matlab files, add it also to Matlab search path
## Instruction for Trnsys Matlab connection using Type 155
* Run the setup file (trnsys-matlab.exe)
* Select your Trnsys installation folder (%TrnsysRoot%) for the destination to unzip files
* Open your TRNSYS workspace in the compiler (%TrnsysRoot%\WorkSpace\CVF66B\CVF66B.dsw)
* Add %TRNSYSRoot%\TrnLib\SEL\SourceCode\Type155.f90 to the "Types" folder under TRNLIB Files (right click, select add, browse to the file)
* Add libeng.lib, libmat.lib and libmx.lib to the "modules" folder under TRNLIB files. Those files are located in %MatlabRoot%\extern\lib\win32\digital\df60\ (%MatlabRoot% is your Matlab installation folder, e.g. C:\Matlab6p5)
* Rebuild Trnlib (F7)
* Launch IISiBat and open %TRNSYSRoot%\TrnLib\SEL\Examples\trnsys-matlab\trnsys-matlab.tpf
* Press F8 to run the simulation. Click on "yes" to exit the online plotter and Matlab will create a new figure
