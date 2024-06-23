package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.LinkedHashMap;
import java.util.Set;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCellStyle;
import org.apache.poi.hssf.usermodel.HSSFFont;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.xssf.usermodel.XSSFCellStyle;
import org.apache.poi.xssf.usermodel.XSSFFont;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.QuantificationThread;
import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.vos.DoubleBondPositionVO;
import at.tugraz.genome.lda.vos.QuantVO;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

public class TargetListFPAnalysis
{
	private static final String IN_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\FP_analysis\\SILDA_II_a_original_FPs.xlsx";
	private static final String OUT_PATH = "D:\\Collaborator_Files\\SILDA\\SILDA_final\\FP_analysis\\out\\SILDA_II_a_original_FPs_analysis\\";
	
	public static void main(String[] args)
  {
  	try
  	{
  		new TargetListFPAnalysis();
  	}
  	catch (Exception ex)
  	{
  		ex.printStackTrace();
  	}
  }
	
	private TargetListFPAnalysis()
	{
		try (XSSFWorkbook inWorkbook = new XSSFWorkbook(IN_PATH);)
		{
			int molSpeciesColumn = -1;
			int commentColumn = -1;
			boolean foundColumns = false;
			int totalIdentifications = 0;
			LinkedHashMap<String,Integer> reasonCount = new LinkedHashMap<String,Integer>();
			
	    for (int sheetNumber = 0; sheetNumber!=inWorkbook.getNumberOfSheets(); sheetNumber++){
	      Sheet sheet = inWorkbook.getSheetAt(sheetNumber);
	      if (sheet.getSheetName().equalsIgnoreCase("SM") || sheet.getSheetName().equalsIgnoreCase("Cer")) continue;
	      for (int rowCount=0;rowCount!=(sheet.getLastRowNum()+1);rowCount++){
	        Row row = sheet.getRow(rowCount);
	        for (int i=0; row!=null && i!=(row.getLastCellNum()+1);i++){
	          Cell cell = row.getCell(i);
	          String contents = "";
	          Double numeric = null;
	          int cellType = -1;
	          if (cell!=null) cellType = cell.getCellType();
	          if (cellType==Cell.CELL_TYPE_STRING){
	            contents = cell.getStringCellValue();
	            try{ 
	              if (contents!=null)numeric = new Double(contents.replaceAll(",", "."));
	            }catch(NumberFormatException nfx){};
	          }else if (cellType==Cell.CELL_TYPE_NUMERIC || cellType==Cell.CELL_TYPE_FORMULA){
	           numeric = cell.getNumericCellValue();
	           contents = String.valueOf(numeric);
	          }
	          if (contents!=null)
	          {
	          	contents = contents.trim();
	          }
	          if (!foundColumns){
	          	if (contents.equalsIgnoreCase("mol. species")){
		          	molSpeciesColumn = i;
	            }
	          	else if (contents.equalsIgnoreCase("Comment")) {
	          		commentColumn = i;
	          	}
	          	if (molSpeciesColumn > -1 && commentColumn > -1)
	          	{
	          		foundColumns = true;
	          	}
	          }
	          else //parse data
	          {
	          	if (i == molSpeciesColumn && contents.length()>1)
	          	{
	          		totalIdentifications++;
	          	}
	          	else if (i == commentColumn && contents.length()>1)
	          	{
	          		if (!reasonCount.containsKey(contents))
	          		{
	          			reasonCount.put(contents, 0);
	          		}
	          		int before = reasonCount.get(contents);
	          		reasonCount.put(contents, ++before);
	          	}
	          }
	        }
	      }
	    }
	    System.out.println("done");
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}
	
	public static XSSFCellStyle getHeaderStyle(XSSFWorkbook wb)
	{
    XSSFCellStyle arial12style = wb.createCellStyle();
    XSSFFont arial12font = wb.createFont();
    arial12font.setBoldweight(HSSFFont.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(HSSFCellStyle.ALIGN_CENTER);
    return arial12style;
  }
}
