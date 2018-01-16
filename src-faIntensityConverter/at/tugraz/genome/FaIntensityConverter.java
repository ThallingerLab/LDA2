package at.tugraz.genome;

import java.io.File;
import java.io.FileNotFoundException;

import at.tugraz.genome.lda.exception.ExcelInputFileException;

public class FaIntensityConverter
{

  public FaIntensityConverter(){
    convertLDAFilesInCurrentDirectory();
  }
  
  
  public static void main(String[] args)
  {
    new FaIntensityConverter();
    
  }

  private void convertLDAFilesInCurrentDirectory(){
    File dir = new File(System.getProperty("user.dir"));
    File[] files = dir.listFiles();
    for (File file : files){
      if (!file.getName().endsWith(".xlsx")) continue;
      LDAToFaConverter converter = new LDAToFaConverter(file.getAbsolutePath());
      try {
        converter.convert();
      } catch (ExcelInputFileException | FileNotFoundException e) {
        System.out.println("---------------------------------------------------");
        System.out.println("There is something wrong with the file: "+file.getName());
        e.printStackTrace();
      }
    }
    
  }
}
