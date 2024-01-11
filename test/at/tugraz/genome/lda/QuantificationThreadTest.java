package at.tugraz.genome.lda;

import static org.junit.jupiter.api.Assertions.*;
import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.AfterEach;

import at.tugraz.genome.lda.exception.SettingsException;
/**
 * 
 * Junit Test class for QuantificationThread.
 * 
 * @author Leonida M. Lamp
 * 
 */
public class QuantificationThreadTest
{
  //The path to the test folder
//  String folderPath = "C:/Users/lample/Collaborator_Files/SILDA/";
//  String folderPath = "E:/Development/LabeledDB/";
//  String folderPath = "C:/Kathi/";
  
  //The sub paths to the test files
//  private static Stream<String> filePaths() {
//    return Stream.of(
//        "20210910_SILDA3/blanks/"
//        "Masslist_Ganglio_positive_all.xlsx"
//      "n-Masslist.xlsx"
//      "massLists/negative.xlsx"
//    );
//  }
  
  QuantificationThread quantThread;
  
  @BeforeEach
  void init()
  {
//    String selectedRawFile = "C:/Users/lample/Collaborator_Files/SILDA/20210910_SILDA3/cells_324_with_A_labeled/200805_324a.chrom";
//    String selectedQuantFile = "C:/Users/lample/Collaborator_Files/SILDA/massLists/negative_A.xlsx";
//    String resultFile = "C:/Users/lample/Collaborator_Files/SILDA/tests/test_200805_324a_negative_A.xlsx";
    String selectedRawFile = "E:/Development/testDB/200803_214a.chrom";
    String selectedQuantFile = "E:/Development/testDB/n-Masslist.xlsx";
//    String selectedQuantFile = "E:/Development/testDB/negative_A.xlsx";
    String resultFile = "E:/Development/testDB/test_200803_214a_n-Masslist_test.xlsx";
    float standardBeforeTolerance = 5f;
    float standardAfterTolerance = 5f;
    int standardAmountOfIsotopes = 2;
    int standardIsotopesMustMatch = 1;
    boolean standardSearchUnknownTime = true;
    float standardBasePeakCutoff = 0.1f;
    float standardRTShift = 0f;
    int numberOfProcessors = 15;
    boolean ionMode = true;
    boolean cli = false;
    
    quantThread = new QuantificationThread(
        selectedRawFile, selectedQuantFile, resultFile,
        standardBeforeTolerance, standardAfterTolerance, 
        standardAmountOfIsotopes, standardIsotopesMustMatch,
        standardSearchUnknownTime,
        standardBasePeakCutoff, standardRTShift, 
        numberOfProcessors, 
        ionMode, cli);
    
    try {
      Settings.saveMachineSettings("OrbiTrap_exactive_SILDA");
      Settings.saveFragmentationSettings("-30_JA_A", "none");
    } catch (SettingsException e) {e.printStackTrace();} 
  }
  
  
  @AfterEach
  void teardown()
  {
    quantThread = null;
  }
  
  @Test
  @DisplayName("Quantifies stuff")
  void runQuantificationThreadTest()
  {
    //when
    int currentLipidCount = 0;
    quantThread.start();
    while (!quantThread.finished()){
      if (quantThread.getCurrentLipidCount()>currentLipidCount){
        System.out.println(String.format("Quantifying %s (%s/%s)", 
            quantThread.getCurrentLipid(), quantThread.getCurrentLipidCount(), quantThread.getTotalAmountOfLipids()));
        currentLipidCount = quantThread.getCurrentLipidCount();
      }
      try { this.wait(1000); } 
      catch ( Exception ex ) {}
    }
    
    //then
    assertTrue(quantThread.finished());
  }
  
  

//  @ParameterizedTest
//  @MethodSource("filePaths")
//  @DisplayName("Parses Quant Files with and without omega-Info.")
//  void parseQuantExcelFileTest(String XLSFile)
//  {
//    //given
//    final String testPath = folderPath+XLSFile;
//    if (new File(testPath).exists() == false) {
//      fail(String.format("The file %s does not exist!", testPath));
//    }
//    
//    //when
//    try
//    {
////      Vector quantValues = QuantificationThread.parseQuantExcelFile(testPath,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f, true);
//      QuantificationThread.parseQuantExcelFile(testPath,  0f, 0f, 0, 0, true, 0f, 0f, 0f, 0f, true);
//    } 
//    catch(Exception ex)
//    {
//      ex.printStackTrace();
//    }  
//  }
}
