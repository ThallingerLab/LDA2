# LDA2

Lipid Data Analyzer 2



Lipid Data Analyzer 2 is released under a GNU GENERAL PUBLIC LICENSE Version 3.

Licensing details can be found in the LICENSE document of this folder.


Project compilation:
For the compilation the following software has to be installed in the following sequence:
Java 8 JDK	http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
Apache Ant	https://ant.apache.org/bindownload.cgi
Java 3D:	https://java3d.java.net/binary-builds.html

After installation and setting the corresponding environment variables vor Java and Ant, the project
can be compiled by opening a command line window, changing the directory to the folder of the
LDA source code (e.g. D:\Development\LipidDataAnalyzer) and typing the following command:
ant

In the toInstall directory of the source code folder, the LipidDataAnalyzer.jar will be created.
The application can be started by executing the LipidDataAnylzer.bat.



Code for MS/MS spectra interpretation:
The code primarily responsible for the MS/MS spectra interpretation can be found in
src/at/tugraz/genome/lda/msn

Scripts for e.g. evaluation of experiment 3 and test cases can be found in:
test/at/tugraz/genome
