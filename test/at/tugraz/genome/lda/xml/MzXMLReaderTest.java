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
import java.io.File;
import java.lang.reflect.Field;

import at.tugraz.genome.lda.swing.Range;
import at.tugraz.genome.maspectras.quantification.CgDefines;

/**
 * 
 * Junit Test class for the MzXMLReader.
 * 
 * @author Leonida M. Lamp
 * 
 */
class MzXMLReaderTest
{
  //To test MzXMLReader class methods NOT initialized by the interface XMLSpectraReader
  MzXMLReader readerMzXML;
 
  //To test methods initialized by the interface XMLSpectraReader
  XMLSpectraReader readerMzXMLInterface;
  
  //The path to the test folder
  String folderPath = "E:/Development/testMzXML/";
  
  //The sub paths to the test files
  private static Stream<String> filePaths() {
    return Stream.of(
//      "absciex/Data20141110_versch_CE_pos+neg1-022b_ref.mzXML",
//      "agilent/Sample1_pos_01_ref.mzXML",
//      "bruker/Mix25_RD1_01_1360_ref.mzXML",
      "thermo/Schev_Meth_test_mono-iso-CID_120731120418_ref.mzXML"
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
    
    readerMzXML = new MzXMLReader(adders, true);
    readerMzXMLInterface = new MzXMLReader(adders, true);    
  }
  
  
  @AfterEach
  void teardown()
  {
    readerMzXML = null;
    readerMzXMLInterface = null;
  }
  

  @ParameterizedTest
  @MethodSource("filePaths")
  @DisplayName("Translates mzXML files without throwing an exception.")
  void translateToChromatogramsTest(String XMLFile)
  {
    //given
    final String testPath = folderPath+XMLFile;
    if (new File(testPath).exists() == false) {
      fail(String.format("The file %s does not exist!", testPath));
    }
    String fileType = "mzXML";
    int maxMBForChromTranslation = 2000;
    int numberOfThreads = 7;
    int multiplicationFactorForInt = 1000;
    int lowestResulution = 1;
    boolean msms = true;
    
    //when
    RawToChromTranslator mzXMLTranslator = new RawToChromTranslator(
        testPath,fileType,maxMBForChromTranslation,numberOfThreads,multiplicationFactorForInt,lowestResulution,msms);
    
    //then   
    try{
      mzXMLTranslator.translateToChromatograms();
    } catch (Exception ex) {
      fail(String.format("An exception occurred while translating the file %s.", testPath));
      ex.printStackTrace();
    }
  }
  
  
  @ParameterizedTest
  @MethodSource("filePaths")
  @DisplayName("Reads mzXML files.")
  void readFileTest(String XMLFile)
  {
    //given
    final String testPath = folderPath+XMLFile;
    if (new File(testPath).exists() == false) {
      fail(String.format("The file %s does not exist!", testPath));
    }
    
    //when
    int[] maximaBefore = {readerMzXMLInterface.getHighestMz(), readerMzXMLInterface.getLowestMz()}; 
    try{
      readerMzXMLInterface.ReadFile(testPath, false);
    } catch (Exception ex) {
      fail(String.format("An exception occurred while translating the file %s.", testPath));
      ex.printStackTrace();
    }
    
    //then
    int[] maximaAfter = {readerMzXMLInterface.getHighestMz(), readerMzXMLInterface.getLowestMz()};
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
    int[] maximaBefore = {readerMzXMLInterface.getHighestMz(), readerMzXMLInterface.getLowestMz()}; 
    try{
      readerMzXMLInterface.ReadFile(testPath, true);
    } catch (Exception ex) {
      fail(String.format("An exception occurred while translating the file %s.", testPath));
      ex.printStackTrace();
    }
    
    //then
    int[] maximaAfter = {readerMzXMLInterface.getHighestMz(), readerMzXMLInterface.getLowestMz()};
    assertFalse(Arrays.equals(maximaBefore, maximaAfter));
  }
  
  
  @Test
  @DisplayName("Gets maxRange_ correctly.")
  void getMaxRangeTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given   
    Range testRange = new Range(1,10);
    final Field field = readerMzXML.getClass().getSuperclass().getDeclaredField("maxRange_");
    field.setAccessible(true);
    field.set(readerMzXML, testRange);

    //when
    final Range result = readerMzXML.getMaxRange();

    //then
    assertEquals(testRange, result);
  }
  
  
  @Test
  @DisplayName("Gets the constant TAG_RUN correctly.")
  void getMsRunNameTest() throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final String expected = "msRun";

    //when
    final String result = readerMzXML.getTagRun();

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
    readerMzXMLInterface.setAdders(newAdder);
    
    //then
    final Field field = readerMzXMLInterface.getClass().getSuperclass().getDeclaredField("adders_");
    field.setAccessible(true);
    assertEquals(newAdder, field.get(readerMzXMLInterface));
  }
  
  
  @ParameterizedTest
  @ValueSource(booleans = {true, false})
  @DisplayName("Gets parseMsMs_ correctly.")
  void getParseMsMsTest(boolean parseMsMs) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzXML.getClass().getSuperclass().getDeclaredField("parseMsMs_");
    field.setAccessible(true);
    field.set(readerMzXML, parseMsMs);

    //when
    final boolean result = readerMzXML.getParseMsMs();

    //then
    assertEquals(parseMsMs, result);
  }

  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Gets multiplicationFactorForInt_ correctly.")
  void getMultiplicationFactorForIntTest(int multiplicationFactorForInt) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzXML.getClass().getSuperclass().getDeclaredField("multiplicationFactorForInt_");
    field.setAccessible(true);
    field.set(readerMzXML, multiplicationFactorForInt);

    //when
    final int result = readerMzXML.getMultiplicationFactorForInt();

    //then
    assertEquals(multiplicationFactorForInt, result);
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Sets highestMz_ correctly.")
  void setHighestMzTest(int highestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //when
    readerMzXML.setHighestMz(highestMz);
    
    //then
    final Field field = readerMzXML.getClass().getSuperclass().getDeclaredField("highestMz_");
    field.setAccessible(true);
    assertEquals(highestMz, field.get(readerMzXML));
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Gets highestMz_ correctly.")
  void getHighestMzTest(int highestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzXMLInterface.getClass().getSuperclass().getDeclaredField("highestMz_");
    field.setAccessible(true);
    field.set(readerMzXMLInterface, highestMz);

    //when
    final int result = readerMzXMLInterface.getHighestMz();

    //then
    assertEquals(highestMz, result);
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Sets lowestMz_ correctly.")
  void setLowestMzTest(int lowestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //when
    readerMzXML.setLowestMz(lowestMz);
    
    //then
    final Field field = readerMzXML.getClass().getSuperclass().getDeclaredField("lowestMz_");
    field.setAccessible(true);
    assertEquals(lowestMz, field.get(readerMzXML));
  }
  
  
  @ParameterizedTest
  @ValueSource(ints = {0, 5, -2021, Integer.MAX_VALUE})
  @DisplayName("Gets lowestMz_ correctly.")
  void getLowestMzTest(int lowestMz) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzXMLInterface.getClass().getSuperclass().getDeclaredField("lowestMz_");
    field.setAccessible(true);
    field.set(readerMzXMLInterface, lowestMz);

    //when
    final int result = readerMzXMLInterface.getLowestMz();

    //then
    assertEquals(lowestMz, result);
  }
  
  
  @ParameterizedTest
  @ValueSource(booleans = {true, false})
  @DisplayName("Sets polaritySwitching_ correctly.")
  void setPolaritySwitchingTest(boolean polaritySwitching) throws NoSuchFieldException, IllegalAccessException
  {
    //when
    readerMzXML.setPolaritySwitching(polaritySwitching);
    
    //then
    final Field field = readerMzXML.getClass().getSuperclass().getDeclaredField("polaritySwitching_");
    field.setAccessible(true);
    assertEquals(polaritySwitching, field.get(readerMzXML));
  }
  
  
  @ParameterizedTest
  @ValueSource(booleans = {true, false})
  @DisplayName("Gets polaritySwitching_ correctly.")
  void usesPolaritySwitchingTest(boolean polaritySwitching) throws NoSuchFieldException, IllegalAccessException
  {
    //given
    final Field field = readerMzXMLInterface.getClass().getSuperclass().getDeclaredField("polaritySwitching_");
    field.setAccessible(true);
    field.set(readerMzXMLInterface, polaritySwitching);

    //when
    final boolean result = readerMzXMLInterface.usesPolaritySwitching();

    //then
    assertEquals(polaritySwitching, result);
  }
  

  
  
}
