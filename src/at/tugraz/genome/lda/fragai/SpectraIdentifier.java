package at.tugraz.genome.lda.fragai;

import java.io.File;
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
import at.tugraz.genome.lda.RawToMzxmlThread;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.lda.xml.AbstractXMLSpectraReader;
import at.tugraz.genome.maspectras.utils.StringUtils;

public class SpectraIdentifier extends Thread
{
	private JTextField rawDirectoryJTextField_;
	
	
	public SpectraIdentifier()
	{
		rawDirectoryJTextField_ = new JTextField(62);
		rawDirectoryJTextField_.setMinimumSize(rawDirectoryJTextField_.getPreferredSize());
		rawDirectoryJTextField_.setToolTipText(TooltipTexts.QUANTITATION_BATCH_RAW_FILE);
		rawDirectoryJTextField_.setText("D:\\Collaborator_Files\\Kathi\\Paper3\\LDA_extension\\Description\\Gangliosides_targets\\TestChrom");
	}
	
	public void identifySpectra()
	{
		translateAllToChrom();
		
		
		
	}
	
	
	
	
	
	
	public void translateAllToChrom()
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
