# LDA2 featuring LC=CL

Lipid Data Analyzer 2 LC=CL extension



Lipid Data Analyzer 2 is released under a GNU GENERAL PUBLIC LICENSE Version 3.

Licensing details can be found in the LICENSE document of this folder.


Project compilation:
For the compilation the following software has to be installed in the following sequence:
Java 8 JDK	http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
Apache Ant	https://ant.apache.org/bindownload.cgi

After installation and setting the corresponding environment variables vor Java and Ant, the project
can be compiled by opening a command line window, changing the directory to the folder of the
LDA source code (e.g. D:\Development\LipidDataAnalyzer) and typing the following command:
ant

In the toInstall directory of the source code folder, the LipidDataAnalyzer.jar will be created.
The application can be started by executing the LipidDataAnylzer.bat.

For starting LDA out of Eclipse, right-click on the LDA project and select "Properties". In the
appearing pop-up, click in the left menu tree on "Java Build Path", select the tab "Libaries",
and click on "Add JARs...". Navigate to the "natives" directory in your project and open the
corresponding directory (linux64, mac or windows64), and all jar files in this directory to your
library path.


Code for LC=CL:
The code primarily responsible for the LC=CL can be found in
src/at/tugraz/genome/lda/target


