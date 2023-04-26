/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2021 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
 * DO NOT ALTER OR REMOVE COPYRIGHT NOTICES OR THIS FILE HEADER. 
 *  
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * by the Free Software Foundation, either version 3 of the License, or 
 * (at your option) any later version.
 *  
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details. 
 *  
 * You should have received a copy of the GNU General Public License 
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 * Please contact lda@genome.tugraz.at if you need additional information or 
 * have any questions.
 */

package at.tugraz.genome.lda.xml;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;
import org.junit.jupiter.params.provider.ValueSource;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.AfterEach;

import java.util.Arrays;
import java.util.Hashtable;
import java.util.stream.Stream;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.lang.reflect.Field;

import at.tugraz.genome.lda.MzxmlToChromThread;
import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.maspectras.quantification.CgDefines;
import at.tugraz.genome.maspectras.utils.StringUtils;

/**
 * 
 * Junit Test class for the MzMLReader.
 * 
 * @author Leonida M. Lamp
 * 
 */
class MzMLReaderTest
{
  //To test MzMLReader class methods NOT initialized by the interface XMLSpectraReader
  MzMLReader readerMzML;
  
  //To test methods initialized by the interface XMLSpectraReader
  XMLSpectraReader readerMzMLInterface;
  
  //The path to the test folder
  String folderPath = "D:\\Collaborator_Files\\testMzXML\\";
  
  //The sub paths to the test files
  private static Stream<String> filePaths() {
    return Stream.of(
      "absciex/testFile.mzML",
      "agilent/testFile.mzML",
      "thermo/testFile.mzML",
      "Kathi_Melbourne/testFile.mzML",
      "bruker/testFile.mzML"
    );
  }
    
  
  @BeforeEach
  void init()
  {
    final Hashtable<Integer,RawToChromThread> translators = new Hashtable<Integer,RawToChromThread>();
    final int numberOfThreads = 3;
    final String[] directories = {"Hello", "World"};
    final int lowerThreshold = 0;
    final int upperThreshold = AbstractXMLSpectraReader.ONE_MILLION;
    
    for (int i=0; i!=numberOfThreads; i++){
      RawToChromThread thread = new RawToChromThread(false,CgDefines.mzMultiplicationFactorForInt,CgDefines.lowestResolution,3);  
      thread.setRequiredInformation(directories, lowerThreshold, upperThreshold);
      translators.put(i,thread);
    }   

    AddScan[] adders = new AddScan[translators.size()];
    for (int i=0;i!= translators.size();i++){
      adders[i] = translators.get(i);
    }
    
    readerMzML = new MzMLReader(adders, true);
    readerMzMLInterface = new MzMLReader(adders, true);    
  }
  
  
  @AfterEach
  void teardown()
  {
    readerMzML = null;
    readerMzMLInterface = null;
  }
  
  
  @ParameterizedTest
  @MethodSource("filePaths")
  @DisplayName("Translates mzML files to chrom as expected.")
  void translateToChromatogramsTest(String XMLFile)
  {
  	//given
  	final String testPath = folderPath+XMLFile;
  	File testFile = new File(testPath);
  	assertTrue(testFile.exists(), String.format("The file %s does not exist!", testPath.toString()));
  	String fileWOExtension = testPath.substring(0,testPath.lastIndexOf("."));
  	File chromFolder = new File(fileWOExtension+".chrom");
  	if (chromFolder.exists()) chromFolder.delete(); //make sure a fully new chrom folder is written.
  	File parent = testFile.getParentFile();
  	File chromFolderReference = new File(parent.toString()+"\\reference.chrom");
  	assertTrue(chromFolderReference.exists(), String.format("The required reference chrom file (%s) does not exist or has a different name!", chromFolderReference.toString()) );
  	String[] subFilesReference = chromFolderReference.list();
  	
    //when
    try{
    	int numProcessors = 1;
//      int numProcessors = Runtime.getRuntime().availableProcessors();
//      numProcessors = numProcessors < 1 ? 1 : numProcessors;
      MzxmlToChromThread thread = new MzXMLToChromThreadTest(testPath,numProcessors);
      thread.run();
      
    } catch (Exception ex) {
      fail(String.format("An exception occurred while translating the file %s: %s", testPath, ex.getMessage()));
    }
    
    //then
    assertTrue(chromFolder.exists(), String.format("The expected chrom file (%s) does not exist or has a different name!", chromFolder.toString()) );
    assertEquals(subFilesReference.length, chromFolder.list().length, String.format("The written chrom file (%s) does not contain the expected number of files.", chromFolder.toString()) );
    
    File headerFile = new File(StringUtils.getChromFilePaths(testPath.toString())[1]);
    assertTrue(headerFile.exists(), "The expected header file does not exist!");
    for (int i=0; i<chromFolder.list().length; i++)
    {
    	String subFile = chromFolder + "\\" + chromFolder.list()[i];
    	String subFileType = subFile.substring(subFile.lastIndexOf(".")+1);
    	if (subFileType.contains("chrom") || subFileType.contains("idx"))
    	{
    		String referenceSubFile = subFile.replace("testFile", "reference");
    		assertTrue(new File(referenceSubFile).exists(), String.format("No reference found for (%s)!", subFile) );
    		try 
        {
        	assertEquals(-1L, filesCompareByByte(subFile, referenceSubFile), String.format("The output file (%s) does not match the reference!", subFile) );
        } 
        catch (IOException ex)
        {
        	fail("An IOException occurred trying to compare the output file to the reference: "+ex.getMessage());
        }
    	}
    }
  }
  
  /**
   * Class to mimick the MzxmlToChromThread without relying on LipidomicsConstants
   * @author Leonida M. Lamp
   */
  private class MzXMLToChromThreadTest extends MzxmlToChromThread
  {
		private MzXMLToChromThreadTest(String filePath, int numberOfThreads)
		{
			super(filePath, numberOfThreads);
		}
		
		@Override
		protected void translateToChrom(String filePath, int numberOfThreads) throws Exception{
	    RawToChromTranslator translator = new RawToChromTranslator(filePath,super.getFileType(filePath),200,
	        numberOfThreads,1000,1,true);
	      translator.translateToChromatograms();
	      super.polaritySwitched_ = translator.isPolaritySwitched();
	  }
  }
  
  private long filesCompareByByte(String path1, String path2) throws IOException 
	{
    try (BufferedInputStream fis1 = new BufferedInputStream(new FileInputStream(new File(path1)));
         BufferedInputStream fis2 = new BufferedInputStream(new FileInputStream(new File(path2)))) 
    {
      int ch = 0;
      long pos = 1;
      while ((ch = fis1.read()) != -1) 
      {
        if (ch != fis2.read()) 
        {
            return pos;
        }
        pos++;
      }
      if (fis2.read() == -1) {
        return -1;
      }
      else {
        return pos;
      }
    }
	}
  
  
  @ParameterizedTest
  @MethodSource("filePaths")
  @DisplayName("Reads mzML file overview.")
  void readFileOverviewTest(String XMLFile)
  {
    //given
    final String testPath = folderPath+XMLFile;
    if (new File(testPath).exists() == false) {
      fail(String.format("The file %s does not exist!", testPath));
    }
    
    //when
    int[] maximaBefore = {readerMzMLInterface.getHighestMz(), readerMzMLInterface.getLowestMz()}; 
    try{
    	readerMzMLInterface.ReadFile(testPath, true);
    } catch (Exception ex) {
      fail(String.format("An exception occurred while reading the file %s: %s", testPath, ex.getMessage()));
    }
    
    //then
    int[] maximaAfter = {readerMzMLInterface.getHighestMz(), readerMzMLInterface.getLowestMz()};
    assertFalse(Arrays.equals(maximaBefore, maximaAfter));
  }
  
  
  @ParameterizedTest
  @MethodSource("filePaths")
  @DisplayName("Reads required info for multi-threading.")
  void readOnlyRequiredInfoForMultiThreadingTest(String XMLFile)
  {
    //given
    final String testPath = folderPath+XMLFile;
    if (new File(testPath).exists() == false) {
      fail(String.format("The file %s does not exist!", testPath));
    }
    
    //when
    int[] maximaBefore = {readerMzMLInterface.getHighestMz(), readerMzMLInterface.getLowestMz()}; 
    try{
      readerMzMLInterface.ReadFile(testPath, true);
    } catch (Exception ex) {
      fail(String.format("An exception occurred while reading the file %s.", testPath));
      ex.printStackTrace();
    }
    
    //then
    int[] maximaAfter = {readerMzMLInterface.getHighestMz(), readerMzMLInterface.getLowestMz()};
    assertFalse(Arrays.equals(maximaBefore, maximaAfter));
  }  
  
  
  @Test
  @DisplayName("Gets maxRange_ correctly.")
  void getMaxRangeTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given   
    Range testRange = new Range(1,10);
    final Field field = readerMzML.getClass().getSuperclass().getDeclaredField("maxRange_");
    field.setAccessible(true);
    field.set(readerMzML, testRange);

    //when
    final Range result = readerMzML.getMaxRange();

    //then
    assertEquals(testRange, result);
  }
  
  
  @Test
  @DisplayName("Gets the constant TAG_RUN correctly.")
  void getMsRunNameTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final String expected = "run";

    //when
    final String result = readerMzML.getTagRun();

    //then
    assertEquals(expected, result);
  }

  
  @Test
  @DisplayName("Sets a new AddScan interface correctly.")
  void setAddersTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    AddScan[] newAdder = new AddScan[1];
    
    //when
    readerMzMLInterface.setAdders(newAdder);
    
    //then
    final Field field = readerMzMLInterface.getClass().getSuperclass().getDeclaredField("adders_");
    field.setAccessible(true);
    assertEquals(newAdder, field.get(readerMzMLInterface));
  }
  
  
  @ParameterizedTest
  @ValueSource(booleans = {true, false})
  @DisplayName("Gets parseMsMs_ correctly.")
  void getParseMsMsTest(boolean parseMsMs) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzML.getClass().getSuperclass().getDeclaredField("parseMsMs_");
    field.setAccessible(true);
    field.set(readerMzML, parseMsMs);

    //when
    final boolean result = readerMzML.getParseMsMs();

    //then
    assertEquals(parseMsMs, result);
  }

  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Gets multiplicationFactorForInt_ correctly.")
  void getMultiplicationFactorForIntTest(int multiplicationFactorForInt) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzML.getClass().getSuperclass().getDeclaredField("multiplicationFactorForInt_");
    field.setAccessible(true);
    field.set(readerMzML, multiplicationFactorForInt);

    //when
    final int result = readerMzML.getMultiplicationFactorForInt();

    //then
    assertEquals(multiplicationFactorForInt, result);
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Sets highestMz_ correctly.")
  void setHighestMzTest(int highestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //when
    readerMzML.setHighestMz(highestMz);
    
    //then
    final Field field = readerMzML.getClass().getSuperclass().getDeclaredField("highestMz_");
    field.setAccessible(true);
    assertEquals(highestMz, field.get(readerMzML));
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Gets highestMz_ correctly.")
  void getHighestMzTest(int highestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzMLInterface.getClass().getSuperclass().getDeclaredField("highestMz_");
    field.setAccessible(true);
    field.set(readerMzMLInterface, highestMz);

    //when
    final int result = readerMzMLInterface.getHighestMz();

    //then
    assertEquals(highestMz, result);
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Sets lowestMz_ correctly.")
  void setLowestMzTest(int lowestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //when
    readerMzML.setLowestMz(lowestMz);
    
    //then
    final Field field = readerMzML.getClass().getSuperclass().getDeclaredField("lowestMz_");
    field.setAccessible(true);
    assertEquals(lowestMz, field.get(readerMzML));
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Gets lowestMz_ correctly.")
  void getLowestMzTest(int lowestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzMLInterface.getClass().getSuperclass().getDeclaredField("lowestMz_");
    field.setAccessible(true);
    field.set(readerMzMLInterface, lowestMz);

    //when
    final int result = readerMzMLInterface.getLowestMz();

    //then
    assertEquals(lowestMz, result);
  }
  
  
  @ParameterizedTest
  @ValueSource(booleans = {true, false})
  @DisplayName("Sets polaritySwitching_ correctly.")
  void setPolaritySwitchingTest(boolean polaritySwitching) throws NoSuchFieldException, IllegalAccessException
  {
    //when
    readerMzML.setPolaritySwitching(polaritySwitching);
    
    //then
    final Field field = readerMzML.getClass().getSuperclass().getDeclaredField("polaritySwitching_");
    field.setAccessible(true);
    assertEquals(polaritySwitching, field.get(readerMzML));
  }
  
  
  @ParameterizedTest
  @ValueSource(booleans = {true, false})
  @DisplayName("Gets polaritySwitching_ correctly.")
  void usesPolaritySwitchingTest(boolean polaritySwitching) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzMLInterface.getClass().getSuperclass().getDeclaredField("polaritySwitching_");
    field.setAccessible(true);
    field.set(readerMzMLInterface, polaritySwitching);

    //when
    final boolean result = readerMzMLInterface.usesPolaritySwitching();

    //then
    assertEquals(polaritySwitching, result);
  }
  

  
  
}
