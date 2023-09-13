package at.tugraz.genome.lda.parser;

import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Properties;
import java.util.Vector;
import java.util.concurrent.ThreadLocalRandom;

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import at.tugraz.genome.lda.LipidomicsConstants;
import at.tugraz.genome.lda.export.QuantificationResultExporter;
import at.tugraz.genome.lda.msn.hydroxy.parser.HydroxyEncoding;
import at.tugraz.genome.lda.quantification.LipidParameterSet;
import at.tugraz.genome.lda.quantification.QuantificationResult;
import at.tugraz.genome.maspectras.quantification.CgAreaStatus;
import at.tugraz.genome.maspectras.quantification.CgException;
import at.tugraz.genome.maspectras.quantification.CgProbe;
import at.tugraz.genome.lda.exception.ExcelInputFileException;
import at.tugraz.genome.lda.exception.ExportException;
import at.tugraz.genome.lda.exception.LipidCombinameEncodingException;
import at.tugraz.genome.lda.exception.SettingsException;

import static org.junit.jupiter.api.Assertions.*;

public class LDAResultReaderIT
{
  private final static String DEFAULT_TEST_DIR = "junit\\testfiles";
  
  
  @Test
  @DisplayName("Writing a QuantificationResult excluding MSn data out and reading it in again leaves the Object unchanged.")
  void readWriteTest() {
    String filePath = DEFAULT_TEST_DIR+"\\quantificationResult.xlsx";
    QuantificationResult quantRes1 = null;
    QuantificationResult quantRes2 = null;
    try {
      quantRes1 = createMockQuantificationResult();
    } catch (LipidCombinameEncodingException ex) {
      fail(ex.getMessage());
    }
    try {
      QuantificationResultExporter.writeResultsToExcel(filePath, quantRes1);
    } catch (ExportException ex) {
      fail(ex.getMessage());
    }
    try {
      quantRes2 = LDAResultReader.readResultFile(filePath, new Hashtable<String,Boolean>());
    } catch (ExcelInputFileException ex) {
      fail(ex.getMessage());
    }

    
    assertTrue(quantRes1.getIdentifications().keySet().equals(quantRes1.getIdentifications().keySet()));
    for (String lipidClass : quantRes1.getIdentifications().keySet()) {
      Vector<LipidParameterSet> lipidParameterSets1 = quantRes1.getIdentifications().get(lipidClass);
      Vector<LipidParameterSet> lipidParameterSets2 = quantRes2.getIdentifications().get(lipidClass);
      assertTrue(lipidParameterSets1.size() == lipidParameterSets2.size());
      for (int i = 0; i < lipidParameterSets1.size(); i++) {
        LipidParameterSet param1 = lipidParameterSets1.get(i);
        LipidParameterSet param2 = lipidParameterSets2.get(i);
        assertTrue(param1.equals(param2));
      }
    }
    assertTrue(quantRes1.getConstants().equals(quantRes2.getConstants()));
  }
  
  
  QuantificationResult createMockQuantificationResult() throws LipidCombinameEncodingException {
    List<String> lipidClasses = new ArrayList<String>();
    lipidClasses.add("ClassOne");
    lipidClasses.add("ClassTwo");
    Hashtable<String,Vector<LipidParameterSet>> identifications = new Hashtable<String,Vector<LipidParameterSet>>();
    for (String lipidClass : lipidClasses) {
      Vector<LipidParameterSet> lipidParameterSets = new Vector<LipidParameterSet>();
      int randomNum = ThreadLocalRandom.current().nextInt(1, 4);
      try {
        lipidParameterSets = createMockLipidParameterSets(randomNum);
      } catch (LipidCombinameEncodingException | CgException ex) {
        throw new LipidCombinameEncodingException(ex);
      }
      identifications.put(lipidClass, lipidParameterSets);
    }
    LipidomicsConstants constants = LipidomicsConstants.getInstance();
    Properties properties = new Properties();
    try {
      constants.setVariables(properties);
    } catch (SettingsException ex) {
      fail(ex.getMessage());
    }
    Map<String,Integer> msLevels = new Hashtable<String,Integer>();
    Hashtable<String,Short> faOhEncondings = new Hashtable<String,Short>();
    Hashtable<String,Short> lcbOhEncondings = new Hashtable<String,Short>();
    HydroxyEncoding faHydroxyEncoding = new HydroxyEncoding(faOhEncondings);
    HydroxyEncoding lcbHydroxyEncoding = new HydroxyEncoding(lcbOhEncondings);
    return new QuantificationResult(identifications, constants, msLevels, faHydroxyEncoding, lcbHydroxyEncoding);
  }
  
  Vector<LipidParameterSet> createMockLipidParameterSets(int num) throws LipidCombinameEncodingException, CgException {
    Vector<LipidParameterSet> mockParams = new Vector<LipidParameterSet>();
    for (int i = 0; i < num; i++) {
      LipidParameterSet mockMS1Set = new LipidParameterSet(829.798461914062f, "48", 1, "NH4", 0.0, "C45 H86 N O8 P", "+N1 +H4",1,0);
      CgProbe probe = new CgProbe(0,1);
      probe.AreaStatus = CgAreaStatus.OK;
      probe.Area = 1.66817472E8f;
      probe.AreaError = 398507968f;
      probe.Background = 4764792f;
      probe.Peak = 989.653991699218f;
      probe.LowerValley = 982.716003417968f;
      probe.UpperValley = 1004.40997314453f;
      probe.Mz = 829.798583984375f;
      probe.LowerMzBand = 0.0130615234375f;
      probe.UpperMzBand = 0.0159912109375f;
      probe.isotopeNumber = 0;
      
      try {
        mockMS1Set.AddProbe(probe);
      } catch (CgException ex) {
        throw new CgException(ex.getMessage());
      }
      mockParams.add(mockMS1Set);
    }
    return mockParams; 
  }
  
  
  
//  void speedTest() {
//    String fileFastExcel = DEFAULT_TEST_DIR+"\\fileFastExcel.xlsx";
//    String fileApachePOI = DEFAULT_TEST_DIR+"\\fileApachePOI.xlsx";
//    final int NB_ROWS = 10000;
//    
//    try {
//      long timeMillis_0 = System.currentTimeMillis();
//      writeFastExcel(fileFastExcel, NB_ROWS);
//      long timeMillis_1 = System.currentTimeMillis();
//      writeApachePOI(fileApachePOI, NB_ROWS);
//      long timeMillis_2 = System.currentTimeMillis();
//      double timeFastExcel = (timeMillis_1-timeMillis_0)/1000.0;
//      double timeApachePOI = (timeMillis_2-timeMillis_1)/1000.0;
//      System.out.println(String.format("Time required by fastExcel: %s \n"
//          + "Time required by Apache POI: %s", timeFastExcel, timeApachePOI));
//    } catch (IOException ex) {
//      System.out.println(ex.getMessage());
//    }
//    
//  }
//  
//  void writeFastExcel(String filePath, int NB_ROWS) throws IOException {
//    try (BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(filePath));) {
//      Workbook wb = new Workbook(out, "Perf", "1.0");
//      Worksheet ws = wb.newWorksheet("Sheet 1");
//      for (int r = 0; r < NB_ROWS; ++r) {
//          ws.value(r, 0, r);
//          ws.value(r, 1, Integer.toString(r % 1000));
//          ws.value(r, 2, r / 87.0);
//          ws.value(r, 3, new Date(1549915044));
//      }
//      ws.range(0, 3, NB_ROWS - 1, 3).style().format("yyyy-mm-dd hh:mm:ss").set();
//      wb.finish();
//    }
//  }
//  
//  void writeApachePOI(String filePath, int NB_ROWS) throws IOException {
//    try (BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(filePath));
//        XSSFWorkbook wb = new XSSFWorkbook();) {
//      Sheet ws = wb.createSheet("Sheet 1");
//      CellStyle dateStyle = wb.createCellStyle();
//      dateStyle.setDataFormat(wb.getCreationHelper().createDataFormat().getFormat("yyyy-mm-dd hh:mm:ss"));
//      for (int r = 0; r < NB_ROWS; ++r) {
//          Row row = ws.createRow(r);
//          row.createCell(0).setCellValue(r);
//          row.createCell(1).setCellValue(Integer.toString(r % 1000));
//          row.createCell(2).setCellValue(r / 87.0);
//          Cell c = row.createCell(3);
//          c.setCellStyle(dateStyle);
//          c.setCellValue(new Date(1549915044));
//      }
//      wb.write(out);
//    }
//  }

}
