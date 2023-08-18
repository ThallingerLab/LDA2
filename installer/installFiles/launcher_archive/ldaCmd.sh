CURRENTDIR=`dirname $0`
cd $CURRENTDIR
java -Xms512m -Xmx1024m -cp LipidDataAnalyzer.jar at.tugraz.genome.lda.LDACmd "$@"