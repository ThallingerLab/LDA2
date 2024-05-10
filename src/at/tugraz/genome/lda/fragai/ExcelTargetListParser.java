package at.tugraz.genome.lda.fragai;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.dhatim.fastexcel.reader.Cell;
import org.dhatim.fastexcel.reader.CellType;
import org.dhatim.fastexcel.reader.ReadableWorkbook;
import org.dhatim.fastexcel.reader.Row;
import org.dhatim.fastexcel.reader.Sheet;

import at.tugraz.genome.lda.utils.StaticUtils;

public class ExcelTargetListParser
{
	private final static String HEADER_CLASS = "Class";
	private final static String HEADER_SPECIES = "Species";
	private final static String HEADER_SUM_FORMULA = "Sum Formula";
	private final static String HEADER_ADDUCTS= "Adducts";
	private final static String HEADER_RETENTION_TIME = "Retention Time";
	private final static String HEADER_TOLERANCE= "Tolerance (+-)";
	private final static String ADDUCT_START= "{";
	private final static String ADDUCT_END= "}";
	private final static String ADDUCT_REGEX= "\\}\\{";
	private final static int HEADER_ROW = 0;
	private File directory_;
	private List<String> headerTitles_;
  private ArrayList<TargetListEntry> targetListEntries_;
  
  /**
	 * Constructor specifying the directory containing the target list in excel format.
	 * @param directory
	 */
	public ExcelTargetListParser(File directory)
	{
		this.directory_ = directory;
		this.targetListEntries_ = new ArrayList<TargetListEntry>();
		
	}
	
  /**
   * Parses the excel file.
   * @throws IOException		if something went wrong parsing the file
   */
  public void parse() throws IOException
  {
  	String supportedFormat = ".xlsx";
    if (!directory_.getAbsolutePath().endsWith(supportedFormat))
    {
    	throw new IOException(String.format("Only the file format '%s' is supported!", supportedFormat));
    }
    try (InputStream is = new FileInputStream(directory_.getAbsolutePath());
        ReadableWorkbook wb = new ReadableWorkbook(is);
        Stream<Sheet> sheets = wb.getSheets();) {
      			sheets.forEach((s) -> { 
            				try { readSheet(s); } 
            				catch (IOException ex) 
            				{
            					headerTitles_ = null; //this indicates that the parsing was not successful
            				} });     
        }
    if (headerTitles_ == null)
    {
    	throw new IOException(String.format("The file '%s' was not parsed successfully. Ensure the required worksheet(s) and headers are named according to the template and entries are in the correct format.\n", directory_.toString()));
    }
  }
  
  /**
	 * Read the excel sheet.
	 * @param sheet					the excel sheet to be read
	 * @throws IOException	if there is something wrong with the excel sheet.
	 */
  private void readSheet(Sheet sheet) throws IOException
	{
		List<Row> rows = null;
    rows = sheet.read();
    Row headerRow = rows.get(HEADER_ROW);
    this.headerTitles_ = readSheetHeaderTitles(headerRow);
    List<Row> contentRows = rows.subList(HEADER_ROW+1, rows.size());
    readContentRows(contentRows);
	}
  
  /**
   * Reads all rows except the header.
   * @param contentRows
   */
  private void readContentRows(List<Row> contentRows) throws IOException
  {
  	if (headerTitles_.contains(HEADER_CLASS) && 
  			headerTitles_.contains(HEADER_SPECIES) && 
  			headerTitles_.contains(HEADER_SUM_FORMULA) &&
  			headerTitles_.contains(HEADER_ADDUCTS) && 
  			headerTitles_.contains(HEADER_RETENTION_TIME) && 
  			headerTitles_.contains(HEADER_TOLERANCE)
  			)
		{
			for (Row row : contentRows) 
			{
	      List<Cell> cells = collectRowCells(row);
	      String lipidClass = null;
	      String species = null;
	    	Hashtable<String,Integer> sumFormula = null;
	    	ArrayList<Adduct> adducts = null;
	    	Double retentionTime = null;
	    	Double tolerance = null;
	    	
	    	int index;
	      String rawValue;
	      
	      for (Cell cell : cells) {
	        index  = cell.getColumnIndex();
	        rawValue = cell.getRawValue();
	        
	        if (index == headerTitles_.indexOf(HEADER_CLASS)) {
	        	lipidClass = rawValue;
	        } else if (index == headerTitles_.indexOf(HEADER_SPECIES)) {
	        	species = rawValue;
	        } else if (index == headerTitles_.indexOf(HEADER_SUM_FORMULA)) {
						try {sumFormula = StaticUtils.categorizeFormula(rawValue);} catch (Exception ex) {}
	        } else if (index == headerTitles_.indexOf(HEADER_ADDUCTS)) {
	        	adducts = new ArrayList<Adduct>();
	        	String[] split = rawValue.split(ADDUCT_REGEX);
	        	for (int i=0;i<split.length;i++)
	        	{
	        		adducts.add(new Adduct(split[i].replace(ADDUCT_START, "").replace(ADDUCT_END, "")));
	        	}
	        } else if (index == headerTitles_.indexOf(HEADER_RETENTION_TIME)) {
	        	retentionTime = Double.parseDouble(rawValue);
	        } else if (index == headerTitles_.indexOf(HEADER_TOLERANCE)) {
	        	tolerance = Double.parseDouble(rawValue);
	        }
	      }
	      if (lipidClass != null && sumFormula != null && adducts != null && retentionTime != null && tolerance != null)
	      {
	      	targetListEntries_.add(new TargetListEntry(lipidClass, species, sumFormula, adducts, retentionTime, tolerance));
	      }
	      else
	      {
	      	throw new IOException(String.format("The content row number %s does not contain all required entries! \n", row.getRowNum()));
	      }
	    }
		}
		else
		{
			throw new IOException("The target list does not contain all required headers.");
		}
  }
  
  /**
   * Parses the header of an Excel sheet
   * @param headerRow 	row number of the header
   * @return List of the header titles
   */
  private List<String> readSheetHeaderTitles(Row headerRow) {
    try (Stream<Cell> cells = headerRow.stream();) {
      return cells.map((c) -> (!(c==null || c.getType().equals(CellType.ERROR)) ? c.getText() : "null")).collect(Collectors.toList());
    } 
  }
  
  /**
   * Collects all cells of a row that are neither null nor contain erroneous values.
   * @param row
   * @return
   */
  private List<Cell> collectRowCells(Row row)
  {
  	return row.stream().filter((c) -> !(c==null || c.getType().equals(CellType.ERROR))).collect(Collectors.toList());
  }
  
  public ArrayList<TargetListEntry> getTargetListEntries()
  {
  	return this.targetListEntries_;
  }
  
  
}