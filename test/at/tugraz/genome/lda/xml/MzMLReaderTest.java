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
  String folderPath = "E:/Development/testMzXML/";
  
  //The sub paths to the test files
  private static Stream<String> filePaths() {
    return Stream.of(
//      "absciex/Data20141110_versch_CE_pos+neg1_proper.mzML",
//      "agilent/Sample1_pos_01_proper.mzML",
//      "bruker/Mix25_RD1_01_1360_proper.mzML",
      "thermo/Schev_Meth_test_mono-iso-CID_120731120418_proper.mzML"
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
  @DisplayName("Translates mzML files without throwing an exception.")
  void translateToChromatogramsTest(String XMLFile)
  {
    //given
    final String testPath = folderPath+XMLFile;
    if (new File(testPath).exists() == false) {
      fail(String.format("The file %s does not exist!", testPath));
    }
    String fileType = "mzML";
    int maxMBForChromTranslation = 2000;
    int numberOfThreads = 7;
    int multiplicationFactorForInt = 1000;
    int lowestResulution = 1;
    boolean msms = true;
    
    //when
    RawToChromTranslator mzMLTranslator = new RawToChromTranslator(
        testPath,fileType,maxMBForChromTranslation,numberOfThreads,multiplicationFactorForInt,lowestResulution,msms);
    
    //then   
    try{
      mzMLTranslator.translateToChromatograms();
    } catch (Exception ex) {
      fail(String.format("An exception occurred while translating the file %s.", testPath));
      ex.printStackTrace();
    }
  }
  
  
  @ParameterizedTest
  @MethodSource("filePaths")
  @DisplayName("Reads mzML files.")
  void readFileTest(String XMLFile)
  {
    //given
    final String testPath = folderPath+XMLFile;
    if (new File(testPath).exists() == false) {
      fail(String.format("The file %s does not exist!", testPath));
    }
    
    //when
    int[] maximaBefore = {readerMzMLInterface.getHighestMz(), readerMzMLInterface.getLowestMz()}; 
    try{
      readerMzMLInterface.ReadFile(testPath, false);
    } catch (Exception ex) {
      fail(String.format("An exception occurred while translating the file %s.", testPath));
      ex.printStackTrace();
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
      fail(String.format("An exception occurred while translating the file %s.", testPath));
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
