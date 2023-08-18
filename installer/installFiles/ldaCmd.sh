cd "${installer:sys.installationDir}"
${installer:javaBinLocation} -Xms${installer:xms}m -Xmx${installer:xmx}m -cp LipidDataAnalyzer.jar at.tugraz.genome.lda.LDACmd "$@"