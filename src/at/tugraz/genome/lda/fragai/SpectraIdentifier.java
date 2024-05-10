package at.tugraz.genome.lda.fragai;

import java.io.File;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.swing.JFrame;
import javax.swing.JTextField;

import at.tugraz.genome.lda.BatchQuantThread;
import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.MzxmlToChromThread;
import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.RawToMzxmlThread;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.exception.QuantificationException;
import at.tugraz.genome.lda.quantification.LipidomicsAnalyzer;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.xml.AbstractXMLSpectraReader;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.maspectras.quantification.ChromatogramReader;
import at.tugraz.genome.maspectras.utils.StringUtils;

public class SpectraIdentifier extends Thread
{
	private JTextField rawDirectoryJTextField_;
	
	
	public SpectraIdentifier()
	{
		rawDirectoryJTextField_ = new JTextField(62);
		rawDirectoryJTextField_.setMinimumSize(rawDirectoryJTextField_.getPreferredSize());
		rawDirectoryJTextField_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_RAW_FILE);
		rawDirectoryJTextField_.setText("D:\\Collaborator_Files\\Kathi\\Paper3\\LDA_extension\\data\\LCMS_STD_Glycolipids_data\\GM1");
	}
	
	public ArrayList<SpectrumContainer> identifySpectra(ExcelTargetListParser parser) throws QuantificationException
	{
		ArrayList<SpectrumContainer> spectra = new ArrayList<SpectrumContainer>();
		translateAllToChrom();
		ArrayList<TargetListEntry> entries = parser.getTargetListEntries();
		
		File rawDir = new File(this.rawDirectoryJTextField_.getText());
    if (rawDir.exists()&&rawDir.isDirectory())
    {
    	File[] chromCandidates = rawDir.listFiles();
    	for (int i=0; i!=chromCandidates.length;i++){
        if (chromCandidates[i].isDirectory() && chromCandidates[i].getAbsolutePath().endsWith("chrom"))
        {
        	String[] chromPaths = StringUtils.getChromFilePaths(chromCandidates[i].getAbsolutePath());
        	ChromatogramReader reader = null;
        	LipidomicsAnalyzer analyzer = null;
        	
        	System.out.println(chromPaths[1]);
        	
        	try
        	{
        		reader = new ChromatogramReader(chromPaths[1], chromPaths[2], chromPaths[3],chromPaths[0],LipidomicsConstants.isSparseData(),LipidomicsConstants.getChromSmoothRange());
          	analyzer = new LipidomicsAnalyzer(chromPaths[1],chromPaths[2],chromPaths[3],chromPaths[0],false); //TODO: one for each thread?
        	}
        	catch (CgException ex)
        	{
        		continue; //move on to the next chrom file
        	}
        	
        	QuantificationThread.setAnalyzerProperties(analyzer);
        	analyzer.setGeneralBasePeakCutoff(0f);
//        	Hashtable<Integer,Float> rtTimes = reader.getRetentionTimesOriginal();
        	int msLevel = 1; //TODO probs edit that at some point
        	float mzTolerance = 0.02f;
        	
        	for (TargetListEntry entry : entries)
        	{
        		Hashtable<String,Integer> formula = entry.getSumFormula();
        		float fullMz = 0f;
      			for (String element : formula.keySet())
      			{
      				fullMz += Settings.getElementParser().getElementDetails(element).getMonoMass()*formula.get(element);
      			}
      			System.out.println(fullMz);
        		ArrayList<Adduct> adducts = entry.getAdducts();
        		for (Adduct adduct : adducts)
        		{
        			float adductMz = fullMz;
        			for (String element : adduct.getAddModifier().keySet())
        			{
        				adductMz += Settings.getElementParser().getElementDetails(element).getMonoMass()*adduct.getAddModifier().get(element);
        			}
        			for (String element : adduct.getRemoveModifier().keySet())
        			{
        				adductMz -= Settings.getElementParser().getElementDetails(element).getMonoMass()*adduct.getRemoveModifier().get(element);
        			}
        			
        			float targetMz = Math.abs(adductMz/adduct.getCharge());
        			
//        			System.out.println(targetMz);
//        			String[] rawLines2D = reader.getRawLines(targetMz-mzTolerance, targetMz+mzTolerance, msLevel);
//        			LipidomicsChromatogram cr = extractChromatogram(targetMz-mzTolerance, targetMz+mzTolerance, rawLines2D, rtTimes, mzTolerance, reader.getMultiplicationFactorForInt_()/reader.getLowestResolution_());
//        			CgProbe probe = analyzer.detectPeakThreeD(cr,LipidomicsAnalyzer.findIndexByTime(new Float(entry.getRetentionTime()*60), cr),false,adduct.getCharge(),msLevel);
        			
        			try
        			{
        				CgProbe probe = analyzer.calculatePeakAtExactTimePosition(new Float(entry.getRetentionTime()*60), targetMz, mzTolerance, mzTolerance, adduct.getCharge(), msLevel);
        				Hashtable<Integer,Vector<String>> spectraRaw = reader.getMsMsSpectra(
          					probe.Mz-LipidomicsConstants.getMs2PrecursorTolerance(), probe.Mz+LipidomicsConstants.getMs2PrecursorTolerance(),probe.LowerValley,probe.UpperValley);
          			
          			@SuppressWarnings("rawtypes")
                Vector<Hashtable> rtNrSpectraAndPrecursor = reader.getRtNrSpectrumHash(spectraRaw);
                @SuppressWarnings("unchecked")
                Hashtable<Integer,String> scanNrSpectrumHash = (Hashtable<Integer,String>)rtNrSpectraAndPrecursor.get(0);
                @SuppressWarnings("unchecked")
                Hashtable<Integer,Vector<Double>> scanNrPrecursorHash = (Hashtable<Integer,Vector<Double>>)rtNrSpectraAndPrecursor.get(1);
                @SuppressWarnings("unchecked")
                Hashtable<Integer,Integer> scanNrLevelHash = (Hashtable<Integer,Integer>)rtNrSpectraAndPrecursor.get(2);
                
                spectra.add(new SpectrumContainer(entry, adduct, probe, chromCandidates[i], scanNrSpectrumHash, scanNrPrecursorHash, scanNrLevelHash));
        			}
              catch (CgException ex)
        			{
              	//do nothing, it is expected that some probes will not have any spectra.
        			}
        		}
        	}
        }
    	}
    }
    return spectra;
	}
	
	

	
	
//	private LipidomicsChromatogram extractChromatogram(float start, float stop, String[] rawLines2D, Hashtable<Integer,Float> rts, float mzTolerance, float resolutionFactor){
//    int startIndex = this.getDataIndex(start, start, resolutionFactor);
//    if (startIndex<1)
//      startIndex = 0;
//    int stopIndex = this.getDataIndex(stop, start, resolutionFactor);
//    if (stopIndex>=rawLines2D.length)
//      stopIndex = rawLines2D.length;
//    LipidomicsChromatogram chrom = new LipidomicsChromatogram(rts.size());
//    chrom.LowerMzBand = mzTolerance;
//    chrom.UpperMzBand = mzTolerance;
//    chrom.Mz = (start+stop)/2;
//    for (int i=0;i!=rts.size();i++){
//      chrom.Value[i][0] = (rts.get(new Integer(i))).floatValue();
//      chrom.Value[i][1] = 0;
//    }
//    if (rawLines2D!=null){
//      if (startIndex<0)
//        startIndex=0;
//      for (int i=startIndex; i<stopIndex && i!=rawLines2D.length; i++){
//        if (rawLines2D[i]!=null&&rawLines2D[i].length()>0){
//          ByteBuffer buffer = ByteBuffer.wrap(Base64.decode(rawLines2D[i]));
//          while (buffer.hasRemaining()){
//            int scanNumber = buffer.getInt();
//            float intensity = buffer.getFloat();
//            chrom.Value[scanNumber][1] += intensity;
//          }
//        }
//      }
//    }
//    chrom.Smooth(LipidomicsConstants.getChromSmoothRange(),
//        LipidomicsConstants.getChromSmoothRepeats());
//    chrom.GetMaximumAndAverage();
//    return chrom;
//  }
//	
//	private int getDataIndex(float mzValue, float start, float resolutionFactor){
//    String toCut = String.valueOf(Calculator.roundFloat((mzValue - start)*resolutionFactor,0));
//    return Integer.parseInt(toCut.substring(0,toCut.indexOf(".")));
//  }
	
	
	private void translateAllToChrom()
	{
    if (rawDirectoryJTextField_.getText()!=null&&rawDirectoryJTextField_.getText().length()>0)
    {
    	File rawDir = new File(this.rawDirectoryJTextField_.getText());
      if (rawDir.exists()&&rawDir.isDirectory())
      {
      	Vector<File> rawFiles = extractRawFiles(rawDir);
      	if (rawFiles.size() > 0)
      	{
      		translateAllRaw(rawFiles);
      	}
      	else{
          if (rawFiles.size()==0){
            new WarningMessage(new JFrame(), "Warning", "In the specified raw directory are no quantifiable files.");
            return;
          }
        }
      	Vector<File> mzXMLFiles = extractMzXMLFiles(rawDir);
      	translateAllMzXML(mzXMLFiles);
      }
      else{
        if (!rawDir.exists()||!rawDir.isDirectory()){
          new WarningMessage(new JFrame(), "Warning", "The raw directory does not exist");
        }
      }
    }
    else {
      if (rawDirectoryJTextField_.getText()==null||rawDirectoryJTextField_.getText().length()<1){
        new WarningMessage(new JFrame(), "Warning", "You must specify a directory containing raw, mzXML, mzML or chrom files.");
      }
    }
	}
	
	
	private void translateAllMzXML(Vector<File> mzXMLFiles)
	{
		for (File file : mzXMLFiles)
  	{
			MzxmlToChromThread mzToChromThread = new MzxmlToChromThread(file.getAbsolutePath(),Runtime.getRuntime().availableProcessors());
			mzToChromThread.start();
			while (!mzToChromThread.finished()) {}
			if (mzToChromThread.getErrorString()==null)
			{
				RawToMzxmlThread.deleteMzXMLFiles(file.getAbsolutePath());
			}
  	}
	}
	
	
	private void translateAllRaw(Vector<File> rawFiles)
	{
		ExecutorService threadpool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
		for (File rawFile : rawFiles)
  	{
  		String suffix = rawFile.getAbsolutePath().substring(rawFile.getAbsolutePath().lastIndexOf("."));
  		if (suffix.equalsIgnoreCase(".RAW")||suffix.equalsIgnoreCase(".d")||suffix.equalsIgnoreCase(".wiff"))
      {
  			RawToMzxmlThread thread = startRawToMzXmlThread(rawFile, suffix);
  			threadpool.execute(thread);
      }
  	}
		threadpool.shutdown();
    try { threadpool.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS); } catch (InterruptedException e) {}
	}
	
	
	private RawToMzxmlThread startRawToMzXmlThread(File rawFile, String suffix)
	{
		if ((rawFile.isFile()&& ((Settings.getReadWPath()!=null&&Settings.getReadWPath().length()>0)||(Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0)))||
        (rawFile.isDirectory() && ((suffix.equalsIgnoreCase(".RAW")&&((Settings.getMassWolfPath()!=null&&Settings.getMassWolfPath().length()>0)||((Settings.getMassPlusPlusPath()!=null&&Settings.getMassPlusPlusPath().length()>0))))
                               ||   (suffix.equalsIgnoreCase(".d") &&Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0)))){
      File headerFile = new File(StringUtils.getChromFilePaths(rawFile.getAbsolutePath())[1]);
      File mzXMLFile = new File(rawFile.getAbsolutePath().substring(0,rawFile.getAbsolutePath().length()-suffix.length())+"."+LipidomicsConstants.getIntermediateFileFormat());
      if (!headerFile.exists()&&!mzXMLFile.exists()){
        boolean isMassPlusPlus = false;
        boolean watersMsConvert = false;
        
        String[] params = new String[3];
        if (rawFile.isFile()){
          if (Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0){
            params = BatchQuantThread.getMsConvertParams(rawFile.getAbsolutePath());
          } else if (Settings.getReadWPath()!=null&&Settings.getReadWPath().length()>0){
            params[0] = Settings.getReadWPath();
            params[1] = rawFile.getAbsolutePath();
            params[2] = "p";
          }  
        }
        if (rawFile.isDirectory()){
          if (suffix.equalsIgnoreCase(".RAW")){
            if (LipidomicsConstants.useMsconvertForWaters()) {
              params =BatchQuantThread.getMsConvertParamsWaters(rawFile.getAbsolutePath());
              watersMsConvert = true;
            } else if (Settings.getMassPlusPlusPath()!=null&&Settings.getMassPlusPlusPath().length()>0){
              params = new String[8];
              params[0] = Settings.getMassPlusPlusPath();
              params[1] = "-in";
              params[2] = rawFile.getAbsolutePath();
              params[3] = "-out";
              params[4] = LipidomicsConstants.getIntermediateFileFormat().toLowerCase(); //Not tested yet for mzML, also unsure whether it has to be lower case..
              params[5] = mzXMLFile.getAbsolutePath();
              params[6] = "-sample";
              params[7] = "0";
              if (LipidomicsConstants.isMS2()) isMassPlusPlus = true;
            }else if (Settings.getMassWolfPath()!=null&&Settings.getMassWolfPath().length()>0){
              params = new String[4];
              params[0] = Settings.getMassWolfPath();
              params[1] = "--"+LipidomicsConstants.getIntermediateFileFormat();
              params[2] = rawFile.getAbsolutePath();
              params[3] = mzXMLFile.getAbsolutePath();  
            }
          }else if(suffix.equalsIgnoreCase(".d")){
            if (Settings.getMsConvertPath()!=null&&Settings.getMsConvertPath().length()>0){
              params =BatchQuantThread.getMsConvertParams(rawFile.getAbsolutePath());
            }                    
          }
        }
        return new RawToMzxmlThread(params,isMassPlusPlus,watersMsConvert);
      }
      else {
      	return null;
      }
    }
		else{
      return null;            
    }
	}
	
	
	private Vector<File> extractRawFiles(File rawDir)
	{
		Vector<File> rawFiles = new Vector<File>();
		File[] rawFileCandidates = rawDir.listFiles();
    Hashtable<String,Vector<File>> avoidDuplication = new Hashtable<String,Vector<File>>();
    boolean mzXMLOrChromPresent = false;
    for (int i=0; i!=rawFileCandidates.length;i++){
      if (rawFileCandidates[i].isFile()){
        String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
        String suffix = fileNameAndSuffix[1];
        String fileName = fileNameAndSuffix[0];
        if (suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML)||
        		suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML)||
        		suffix.equalsIgnoreCase("raw")||
        		suffix.equalsIgnoreCase("chrom")||
        		suffix.equalsIgnoreCase("wiff")){
          if (suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML)||suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML)||suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
          Vector<File> theFiles = new Vector<File>();
          if (avoidDuplication.containsKey(fileName)){
            theFiles = avoidDuplication.get(fileName);
          }
          theFiles.add(rawFileCandidates[i]);
          avoidDuplication.put(fileName, theFiles);
        }
      }
      if (rawFileCandidates[i].isDirectory()){
        String[] fileNameAndSuffix = StaticUtils.extractFileNameAndSuffix(rawFileCandidates[i].getAbsolutePath()); 
        String suffix = fileNameAndSuffix[1];
        String fileName = fileNameAndSuffix[0];
        if (suffix.equalsIgnoreCase("raw")|| suffix.equalsIgnoreCase("d") ||suffix.equalsIgnoreCase("chrom")){
          if (suffix.equalsIgnoreCase("chrom")) mzXMLOrChromPresent = true;
          Vector<File> theFiles = new Vector<File>();
          if (avoidDuplication.containsKey(fileName)){
            theFiles = avoidDuplication.get(fileName);
          }
          theFiles.add(rawFileCandidates[i]);
          avoidDuplication.put(fileName, theFiles);
        }
      }
    }
    for (String key : avoidDuplication.keySet()){
      Vector<File> theFiles = avoidDuplication.get(key);
      if (theFiles.size()==1){
        String suffix = StaticUtils.extractFileNameAndSuffix(theFiles.get(0).getAbsolutePath())[1];
        if (!mzXMLOrChromPresent || !suffix.equalsIgnoreCase("wiff"))
          rawFiles.add(theFiles.get(0));
      }else{
        int selectedIndex = -1;
        for (int i=0; i!=theFiles.size();i++){
          File file = theFiles.get(i);
          String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
          if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
          if (suffix.equalsIgnoreCase("chrom")){
            selectedIndex = i;
          }
        }
        if (selectedIndex>-1){
          rawFiles.add(theFiles.get(selectedIndex));
        }else{
          for (int i=0; i!=theFiles.size();i++){
            File file = theFiles.get(i);
            String suffix = file.getAbsolutePath().substring(file.getAbsolutePath().lastIndexOf(".")+1);
            if (mzXMLOrChromPresent && suffix.equalsIgnoreCase("wiff")) continue;
            if (suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML) || suffix.equalsIgnoreCase(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML)){
              rawFiles.add(theFiles.get(i));
            }
          }  
        }
      }
    }
    return rawFiles;
	}
	
	private Vector<File> extractMzXMLFiles(File rawDir)
	{
		Vector<File> files = new Vector<File>();
		File[] mzMLFileCandidates = rawDir.listFiles();
		for (int i=0; i<mzMLFileCandidates.length;i++){
			if (mzMLFileCandidates[i].isFile()){
        if (mzMLFileCandidates[i].getAbsolutePath().endsWith(AbstractXMLSpectraReader.FILE_TYPE_MZ_XML)||mzMLFileCandidates[i].getAbsolutePath().endsWith(AbstractXMLSpectraReader.FILE_TYPE_MZ_ML))
        {
        	files.add(mzMLFileCandidates[i]);
        }
			}
		}
		return files;
	}
	
	
}
