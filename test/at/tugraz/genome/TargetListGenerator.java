package at.tugraz.genome;

import java.io.BufferedOutputStream;
import java.io.FileOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;

import org.apache.poi.hssf.usermodel.HSSFCell;
import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.xssf.usermodel.XSSFRow;
import org.apache.poi.xssf.usermodel.XSSFSheet;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;

import at.tugraz.genome.lda.exception.ChemicalFormulaException;
import at.tugraz.genome.lda.msn.parser.FALibParser;
import at.tugraz.genome.lda.msn.vos.FattyAcidVO;
import at.tugraz.genome.lda.utils.StaticUtils;
import at.tugraz.genome.maspectras.parser.spectrummill.ElementConfigParser;

public class TargetListGenerator
{
	String chainListLabelCombi_ = "D";
	ArrayList<String> lipidClasses_ = new ArrayList<String>(Arrays.asList("PC","O-PC","P-PC","LPC","PE","P-PE","LPE","PS","LPS","PI","LPI","PG","LPG","SM","Cer"));
//	ArrayList<String> lipidClassesGanglio_ = new ArrayList<String>(Arrays.asList("GM1", "GM2", "GM3", "GM4", "GD1", "GD2", "GD3", "GT1","GT2","GT3", "GQ1","GP1", "GH1", "GS1", "GO1", "LacCer"));
	ArrayList<String> lipidClassesGanglio_ = new ArrayList<String>(Arrays.asList(
			"GM1", "GM1-Ac", "GM1-Fuc", "GM2", "GM3", "GM4", 
			"GD1", "GD1-Ac", "GD1-Fuc", "GD1-Fuc-Ac", "GD2", "GD3", "GD3-Ac", 
			"GT1", "GT1-Ac", "GT1-Fuc", "GT2", "GT3", "GT3-Ac", "GT3-Fuc", 
			"GQ1","GP1", "GH1", "GS1", "GS1-Fuc", "GO1", "GQ1-Ac", "GQ1-Fuc", "LacCer"));
	String chainListPath_ = "D:\\Collaborator_Files\\Kathi\\Paper2\\newTL\\Ganglio_chains\\fattyAcidChains.xlsx";
	String lcbListPath_ = "D:\\Collaborator_Files\\Kathi\\Paper2\\newTL\\Ganglio_chains\\dLCB_Ganglio.xlsx";
	String outFile_ = "D:\\Collaborator_Files\\Kathi\\Paper2\\newTL\\Ganglio_negative_FucAc.xlsx";
	
	
//	String chainListPath_ = String.format("D:\\Collaborator_Files\\SILDA\\SILDA_final\\masslists_all_labels\\new_TL_D\\chainLists\\fattyAcidChains_%s.xlsx", chainListLabelCombi_);
//	String outFile_ = String.format("D:\\Collaborator_Files\\SILDA\\SILDA_final\\masslists_all_labels\\new_TL_D\\negative_%s.xlsx", chainListLabelCombi_);
	
	public static void main(String[] args)
  {
		new TargetListGenerator();
  }
	
	private TargetListGenerator()
	{
		generateMassListGanglio(lipidClassesGanglio_, chainListPath_, lcbListPath_, outFile_);
//		generateIsoLabeledMassListNew(lipidClasses_, chainListPath_, outFile_);
	}
	
	
	private void generateMassListGanglio(ArrayList<String> lipidClasses, String chainListPath, String lcbListPath, String outFile)
	{
		try {
			ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFile));
      XSSFWorkbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      for (String lipidClass : lipidClasses)
      {
      	XSSFSheet resultSheet = resultWorkbook.createSheet(lipidClass);
    		int rowCount = 1;
    		XSSFRow headerRow = resultSheet.createRow(rowCount);
        Cell cell = headerRow.createCell(10,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("OH-Number: 2");
        cell = headerRow.createCell(11,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("OH-Range: 2-3");
        cell = headerRow.createCell(12,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("adductInsensitiveRtFilter");
        
      	rowCount = 4;
        headerRow = resultSheet.createRow(rowCount);
        cell = headerRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("Name");
        cell.setCellStyle(headerStyle);
        cell = headerRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("dbs");
        cell.setCellStyle(headerStyle);
        cell = headerRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("C");
        cell.setCellStyle(headerStyle);
        cell = headerRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("H");
        cell.setCellStyle(headerStyle);
        cell = headerRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("O");
        cell.setCellStyle(headerStyle);
        cell = headerRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("N");
        cell.setCellStyle(headerStyle);
        cell = headerRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("M");
        cell.setCellStyle(headerStyle);
        
        final String formH1 = "mass(form[-H] name[-H])";
        final String formH2 = "mass(form[-H2] name[-H2] charge=2)";
        final String formH3 = "mass(form[-H3] name[-H3] charge=3)";
        final String formH4 = "mass(form[-H4] name[-H4] charge=4)";
        final String formH5 = "mass(form[-H5] name[-H5] charge=5)";
        final String formH6 = "mass(form[-H6] name[-H6] charge=6)";
        final String formH7 = "mass(form[-H7] name[-H7] charge=7)";
        final String formH8 = "mass(form[-H8] name[-H8] charge=8)";
        
        int lastColumn = 7;
        if (lipidClass.equals("LacCer"))
        {
        	cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH1);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH2);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GM1") || lipidClass.contains("GM2") || lipidClass.contains("GM3") || lipidClass.contains("GM4"))
        {
        	cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH1);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH2);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GD1") || lipidClass.contains("GD2") || lipidClass.contains("GD3"))
        {
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH2);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GT1") || lipidClass.contains("GT2") || lipidClass.contains("GT3"))
        {
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH2);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH3);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GQ1"))
        {
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH3);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH4);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GP1"))
        {
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH4);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH5);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GS1"))
        {
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH4);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH5);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH6);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GH1"))
        {
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH4);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH5);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH6);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH7);
          cell.setCellStyle(headerStyle);
        }
        else if (lipidClass.contains("GO1"))
        {
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH4);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH5);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH6);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH7);
          cell.setCellStyle(headerStyle);
          cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue(formH8);
          cell.setCellStyle(headerStyle);
        }
        
        cell = headerRow.createCell(++lastColumn,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("tR (min)");
        cell.setCellStyle(headerStyle);
        rowCount++;
        
        int minFattyAcidC = 30;
        int maxFattyAcidC = 51;
        int maxDB=6;
        Hashtable<Integer,Integer> cDbsCombi = new Hashtable<Integer,Integer>();
        
      	for (int i=minFattyAcidC; i<=maxFattyAcidC; i++) 
      	{
          for (int j=0; j<=maxDB; j++) 
          {
            if (j>4 && i<30)
              continue;
            else if (j>6 && i<36)
              continue;
            else if (j>8 && i<38)
              continue;
            else if (j>10 && i<40)
              continue;
            cDbsCombi.put(i, j);
          }
        }
        
        for (int fattyAcidC=minFattyAcidC; fattyAcidC<maxFattyAcidC; fattyAcidC+=1)
        {
          int dbs = 0;
          int dbsMax = cDbsCombi.get(fattyAcidC)+1;
          while (dbs<dbsMax){        
            Hashtable<String,Hashtable<String,Integer>> allElements = new Hashtable<String,Hashtable<String,Integer>>();
            allElements = fillUnlabeledElements(allElements, fattyAcidC, dbs, lipidClass);
              
            for (String current : allElements.keySet())
            {
            	Hashtable<String,Integer> elements = allElements.get(current);
              double massNeutral = 
              		parser.getElementDetails("C").getMonoMass()*elements.get("C")+
              		parser.getElementDetails("Cc").getMonoMass()*elements.get("Cc")+
              		parser.getElementDetails("H").getMonoMass()*elements.get("H")+
              		parser.getElementDetails("D").getMonoMass()*elements.get("D")+
                  parser.getElementDetails("O").getMonoMass()*elements.get("O")+
                  parser.getElementDetails("P").getMonoMass()*elements.get("P")+
                  parser.getElementDetails("N").getMonoMass()*elements.get("N");
              
              XSSFRow outRow = resultSheet.createRow(rowCount);
              cell = outRow.createCell(0,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(current+fattyAcidC);
              cell = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
              cell.setCellValue(":");       
              cell = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(dbs);
              cell = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("C"));
              cell = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("H"));
              cell = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("O"));
              cell = outRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("N"));
              cell = outRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(massNeutral);
              
              for (int i=8;i<=lastColumn;i++)
              {
              	double massAdduct = massNeutral;
              	if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH1))
              	{
              		massAdduct = massNeutral-parser.getElementDetails("h").getMonoMass();
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
              	else if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH2))
              	{
              		massAdduct = (massNeutral-2*parser.getElementDetails("h").getMonoMass())/2;
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
              	else if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH3))
              	{
              		massAdduct = (massNeutral-3*parser.getElementDetails("h").getMonoMass())/3;
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
              	else if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH4))
              	{
              		massAdduct = (massNeutral-4*parser.getElementDetails("h").getMonoMass())/4;
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
              	else if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH5))
              	{
              		massAdduct = (massNeutral-5*parser.getElementDetails("h").getMonoMass())/5;
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
              	else if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH6))
              	{
              		massAdduct = (massNeutral-6*parser.getElementDetails("h").getMonoMass())/6;
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
              	else if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH7))
              	{
              		massAdduct = (massNeutral-7*parser.getElementDetails("h").getMonoMass())/7;
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
              	else if (headerRow.getCell(i).getStringCellValue().equalsIgnoreCase(formH8))
              	{
              		massAdduct = (massNeutral-8*parser.getElementDetails("h").getMonoMass())/8;
              		cell = outRow.createCell(i,HSSFCell.CELL_TYPE_NUMERIC);
                  cell.setCellValue(massAdduct);
              	}
                
              }
              
              rowCount++;
            }
            dbs++;
          }
        }
      }
      resultWorkbook.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
    System.out.println("workbook written!");
	}
	
	private Vector<FattyAcidVO> getChains(String path) throws Exception
	{
		Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
		FALibParser faParser = new FALibParser(path);
    faParser.parseFile();
    Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> fas = faParser.getFattyAcids().get("n");
    for (Integer cAtoms : fas.keySet()) {
      for (Integer dbs : fas.get(cAtoms).keySet()) {
        for (String pref : fas.get(cAtoms).get(dbs).keySet()) {
      	  for (String oxState : fas.get(cAtoms).get(dbs).get(pref).keySet()) {
      	  	chains.add(fas.get(cAtoms).get(dbs).get(pref).get(oxState));
      	  }
        }
      }
    }
    return chains;
	}
	
	
	private void generateIsoLabeledMassListNew(ArrayList<String> lipidClasses, String chainListPath, String outFile) {
    try {
    	ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      FALibParser faParser = new FALibParser(chainListPath);
      faParser.parseFile();
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> fas = faParser.getFattyAcids().get("n");
      Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
      for (Integer cAtoms : fas.keySet()) {
        for (Integer dbs : fas.get(cAtoms).keySet()) {
          for (String pref : fas.get(cAtoms).get(dbs).keySet()) {
        	  for (String oxState : fas.get(cAtoms).get(dbs).get(pref).keySet()) {
        		  chains.add(fas.get(cAtoms).get(dbs).get(pref).get(oxState));
        	  }
          }
        }
      }
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFile));
      XSSFWorkbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      for (String lipidClass : lipidClasses)
      {
      	XSSFSheet resultSheet = resultWorkbook.createSheet(lipidClass);
      	if (lipidClass.equals("SM") || lipidClass.equals("Cer"))
      	{
      		int rowCount = 1;
      		XSSFRow outRow = resultSheet.createRow(rowCount);
          Cell cell = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);
          cell.setCellValue("OH-Number: 2");
          cell = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);
          cell.setCellValue("OH-Range: 2-3");
          cell = outRow.createCell(12,HSSFCell.CELL_TYPE_STRING);
          cell.setCellValue("adductInsensitiveRtFilter");
      	}
      	int rowCount = 4;
        XSSFRow outRow = resultSheet.createRow(rowCount);
        Cell cell = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("Name");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("dbs");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("C");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("Cc");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("H");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("D");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("O");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(8,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("P");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("N");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("M");
        cell.setCellStyle(headerStyle);
        
        if (lipidClass.contains("PC") || lipidClass.equals("SM") || lipidClass.equals("Cer"))
        {
        	cell = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue("mass(form[+HCOO] name[HCOO])");
          cell.setCellStyle(headerStyle);
        }
        else
        {
        	cell = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue("mass(form[-H] name[-H])");
          cell.setCellStyle(headerStyle);
        }
        
        cell = outRow.createCell(12,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("tR (min)");
        cell.setCellStyle(headerStyle);
        rowCount++;
        
        int minFattyAcidC = Integer.MAX_VALUE;
        int maxFattyAcidC = Integer.MIN_VALUE;
        int numFA = 0;
        int maxDB;
        Hashtable<Integer,Integer> cDbsCombi = new Hashtable<Integer,Integer>();
        if (lipidClass.startsWith("L"))
        {
        	minFattyAcidC = 12;
        	maxFattyAcidC = 27;
        	numFA=1;
        	maxDB=8;
        	for (int i=minFattyAcidC; i<=maxFattyAcidC; i++) 
        	{
            for (int j=0; j<=maxDB; j++) 
            {
              if (j>4 && i<16)
                continue;
              else if (j>6 && i<20)
                continue;
              cDbsCombi.put(i, j);
            }
          }
        }
        else
        {
        	minFattyAcidC = 20;
        	maxFattyAcidC = 49;
        	numFA=2;
        	maxDB=12;
        	if (lipidClass.equals("SM") || lipidClass.equals("Cer"))
        	{
        		maxDB=6;
        	}
        	for (int i=minFattyAcidC; i<=maxFattyAcidC; i++) 
        	{
            for (int j=0; j<=maxDB; j++) 
            {
              if (j>4 && i<30)
                continue;
              else if (j>6 && i<36)
                continue;
              else if (j>8 && i<38)
                continue;
              else if (j>10 && i<40)
                continue;
              cDbsCombi.put(i, j);
            }
          }
        }
        
        for (int fattyAcidC=minFattyAcidC; fattyAcidC<maxFattyAcidC; fattyAcidC++)
        {
          int dbs = 0;
          int dbsMax = cDbsCombi.get(fattyAcidC)+1;
          while (dbs<dbsMax){
            Vector<Vector<FattyAcidVO>> filtered = getPossCombis(numFA,fattyAcidC,dbs,chains,new Vector<FattyAcidVO>());          
            Hashtable<String,Hashtable<String,Integer>> labelElements = new Hashtable<String,Hashtable<String,Integer>>();
            labelElements = fillLabelElements(labelElements,fattyAcidC,dbs,lipidClass,filtered);
              
            for (String label : labelElements.keySet())
            {
            	if ((lipidClass.equals("SM") || lipidClass.equals("Cer")) && label.length()>1)
            	{
            		continue;
            	}
            	if (label.length()<1)
            	{
            		continue;
            	}
            	Hashtable<String,Integer> elements = labelElements.get(label);
              double massNeutral = 
              		parser.getElementDetails("C").getMonoMass()*elements.get("C")+
              		parser.getElementDetails("Cc").getMonoMass()*elements.get("Cc")+
              		parser.getElementDetails("H").getMonoMass()*elements.get("H")+
              		parser.getElementDetails("D").getMonoMass()*elements.get("D")+
                  parser.getElementDetails("O").getMonoMass()*elements.get("O")+
                  parser.getElementDetails("P").getMonoMass()*elements.get("P")+
                  parser.getElementDetails("N").getMonoMass()*elements.get("N");
              
              double massAdduct = massNeutral;
              if (lipidClass.contains("PC") || lipidClass.equals("SM") || lipidClass.equals("Cer")) //formate
              {
              	massAdduct = massNeutral-
              			parser.getElementDetails("h").getMonoMass()+
              			parser.getElementDetails("H").getMonoMass()*2+
              			parser.getElementDetails("C").getMonoMass()+
              			parser.getElementDetails("O").getMonoMass()*2;
              }
              else //deprotonated
              {
              	massAdduct = massNeutral-
              			parser.getElementDetails("h").getMonoMass();
              }
              
              outRow = resultSheet.createRow(rowCount);
              cell = outRow.createCell(0,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(label+fattyAcidC);
              cell = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
              cell.setCellValue(":");       
              cell = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(dbs);
              cell = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("C"));
              cell = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("Cc"));
              cell = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("H"));
              cell = outRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("D"));
              cell = outRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("O"));
              cell = outRow.createCell(8,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("P"));
              cell = outRow.createCell(9,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("N"));
              cell = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(massNeutral);
              cell = outRow.createCell(11,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(massAdduct);
              rowCount++;
            }
            dbs++;
          }
        }
      }
      resultWorkbook.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
    System.out.println("workbook written!");
  }
	
	
	private void generateIsoLabeledMassListGosia() {
  	String chainListLabelCombi = "d27_Myristic";
  	ArrayList<String> lipidClasses = new ArrayList<String>(Arrays.asList("PC","PE","P-PE","PS","PI","PG"));
  	String chainListPath = String.format("D:\\Collaborator_Files\\Gosia\\fattyAcidChains_%s.xlsx", chainListLabelCombi);
  	String outFile = String.format("D:\\Collaborator_Files\\Gosia\\negative_%s.xlsx", chainListLabelCombi);
    try {
    	ElementConfigParser parser = new ElementConfigParser("elementconfig.xml");
      parser.parse();
      FALibParser faParser = new FALibParser(chainListPath);
      faParser.parseFile();
      Hashtable<Integer,Hashtable<Integer,Hashtable<String,Hashtable<String,FattyAcidVO>>>> fas = faParser.getFattyAcids().get("n");
      Vector<FattyAcidVO> chains = new Vector<FattyAcidVO>();
      for (Integer cAtoms : fas.keySet()) {
        for (Integer dbs : fas.get(cAtoms).keySet()) {
          for (String pref : fas.get(cAtoms).get(dbs).keySet()) {
        	  for (String oxState : fas.get(cAtoms).get(dbs).get(pref).keySet()) {
        		  chains.add(fas.get(cAtoms).get(dbs).get(pref).get(oxState));
        	  }
            
          }
        }
      }
      BufferedOutputStream out = new BufferedOutputStream(new FileOutputStream(outFile));
      XSSFWorkbook resultWorkbook = new XSSFWorkbook();
      CellStyle headerStyle = getHeaderStyle(resultWorkbook);
      for (String lipidClass : lipidClasses)
      {
      	XSSFSheet resultSheet = resultWorkbook.createSheet(lipidClass);
      	int rowCount = 4;
        XSSFRow outRow = resultSheet.createRow(rowCount);
        Cell cell = outRow.createCell(0,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("Name");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(2,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("dbs");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(3,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("C");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(4,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("Cc");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(5,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("H");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(6,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("D");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(7,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("O");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(8,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("P");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(9,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("N");
        cell.setCellStyle(headerStyle);
        cell = outRow.createCell(10,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("M");
        cell.setCellStyle(headerStyle);
        
        if (lipidClass.contains("PC") || lipidClass.equals("SM"))
        {
        	cell = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue("mass(form[+HCOO] name[HCOO])");
          cell.setCellStyle(headerStyle);
        }
        else
        {
        	cell = outRow.createCell(11,HSSFCell.CELL_TYPE_STRING);      
          cell.setCellValue("mass(form[-H] name[-H])");
          cell.setCellStyle(headerStyle);
        }
        
        cell = outRow.createCell(12,HSSFCell.CELL_TYPE_STRING);
        cell.setCellValue("tR (min)");
        cell.setCellStyle(headerStyle);
        rowCount++;
        
        int minFattyAcidC = Integer.MAX_VALUE;
        int maxFattyAcidC = Integer.MIN_VALUE;
        int numFA = 0;
        Hashtable<Integer,Integer> cDbsCombi = new Hashtable<Integer,Integer>();
        if (lipidClass.equals("LPC") || lipidClass.equals("LPE") || lipidClass.equals("LPS") || lipidClass.equals("SM"))
        {
        	minFattyAcidC = 12;
        	maxFattyAcidC = 27;
        	numFA=1;
        	for (int i=minFattyAcidC; i!=maxFattyAcidC; i++) 
        	{
            for (int j=0; j!=9; j++) 
            {
              if (j>4 && i<16)
                continue;
              else if (j>6 && i<20)
                continue;
              cDbsCombi.put(i, j);
            }
          }
        }
        else
        {
        	minFattyAcidC = 20;
        	maxFattyAcidC = 49;
        	numFA=2;
        	for (int i=minFattyAcidC; i!=maxFattyAcidC; i++) 
        	{
            for (int j=0; j!=13; j++) 
            {
              if (j>4 && i<30)
                continue;
              else if (j>6 && i<36)
                continue;
              else if (j>8 && i<38)
                continue;
              else if (j>10 && i<40)
                continue;
              cDbsCombi.put(i, j);
            }
          }
        }
        
        for (int fattyAcidC=minFattyAcidC; fattyAcidC<maxFattyAcidC; fattyAcidC+=1)
        {
          int dbs = 0;
          int dbsMax = cDbsCombi.get(fattyAcidC)+1;
          while (dbs<dbsMax){
            Vector<Vector<FattyAcidVO>> filtered = getPossCombis(numFA,fattyAcidC,dbs,chains,new Vector<FattyAcidVO>());          
            Hashtable<String,Hashtable<String,Integer>> labelElements = new Hashtable<String,Hashtable<String,Integer>>();
            labelElements = fillLabelElements(labelElements,fattyAcidC,dbs,lipidClass,filtered);
              
            for (String label : labelElements.keySet())
            { 
            	Hashtable<String,Integer> elements = labelElements.get(label);
              double massNeutral = 
              		parser.getElementDetails("C").getMonoMass()*elements.get("C")+
              		parser.getElementDetails("Cc").getMonoMass()*elements.get("Cc")+
              		parser.getElementDetails("H").getMonoMass()*elements.get("H")+
              		parser.getElementDetails("D").getMonoMass()*elements.get("D")+
                  parser.getElementDetails("O").getMonoMass()*elements.get("O")+
                  parser.getElementDetails("P").getMonoMass()*elements.get("P")+
                  parser.getElementDetails("N").getMonoMass()*elements.get("N");
              
              double massAdduct = massNeutral;
              if (lipidClass.contains("PC") || lipidClass.equals("SM")) //formate
              {
              	massAdduct = massNeutral-
              			parser.getElementDetails("h").getMonoMass()+
              			parser.getElementDetails("H").getMonoMass()*2+
              			parser.getElementDetails("C").getMonoMass()+
              			parser.getElementDetails("O").getMonoMass()*2;
              }
              else //deprotonated
              {
              	massAdduct = massNeutral-
              			parser.getElementDetails("h").getMonoMass();
              }
              
              outRow = resultSheet.createRow(rowCount);
              cell = outRow.createCell(0,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(label+fattyAcidC);
              cell = outRow.createCell(1,HSSFCell.CELL_TYPE_STRING);
              cell.setCellValue(":");       
              cell = outRow.createCell(2,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(dbs);
              cell = outRow.createCell(3,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("C"));
              cell = outRow.createCell(4,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("Cc"));
              cell = outRow.createCell(5,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("H"));
              cell = outRow.createCell(6,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("D"));
              cell = outRow.createCell(7,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("O"));
              cell = outRow.createCell(8,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("P"));
              cell = outRow.createCell(9,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(elements.get("N"));
              cell = outRow.createCell(10,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(massNeutral);
              cell = outRow.createCell(11,HSSFCell.CELL_TYPE_NUMERIC);
              cell.setCellValue(massAdduct);
              rowCount++;
            }
            dbs++;
          }
        }
      }
      resultWorkbook.write(out);
      out.close();
    } catch (Exception ex){
      ex.printStackTrace();
    }
    System.out.println("workbook written!");
  }
	
	
	private CellStyle getHeaderStyle(Workbook wb){
    CellStyle arial12style = wb.createCellStyle();
    Font arial12font = wb.createFont();
    arial12font.setBoldweight(Font.BOLDWEIGHT_BOLD);
    arial12font.setFontName("Arial");
    arial12font.setFontHeightInPoints((short)12);
    arial12style.setFont(arial12font);
    arial12style.setAlignment(CellStyle.ALIGN_CENTER);
    return arial12style;
  }
	
	private Hashtable<String,Hashtable<String,Integer>> fillLabelElements(Hashtable<String,Hashtable<String,Integer>>labelElements, 
  		int fattyAcidC, int dbs, String lipidClass, Vector<Vector<FattyAcidVO>> filtered) throws ChemicalFormulaException
  {
  	labelElements = fillUnlabeledElements(labelElements,fattyAcidC,dbs,lipidClass);
  	for (Vector<FattyAcidVO> combs : filtered) 
    {
    	Vector<String> labels = new Vector<String>();
    	Hashtable<String,Integer> classElements = new Hashtable<String,Integer>(labelElements.get(""));
    	for (FattyAcidVO fa : combs) 
    	{
    		if (fa.getPrefix().length()>0) 
    		{
    			labels.add(fa.getPrefix());
    			Hashtable<String,Integer> formula = StaticUtils.categorizeFormula(fa.getFormula(), true);
    			if (formula.containsKey("D"))
    			{
    				classElements.put("D", classElements.get("D")+formula.get("D"));
    				classElements.put("H", classElements.get("H")-formula.get("D"));
    			}
    			if (formula.containsKey("Cc"))
    			{
    				classElements.put("Cc", classElements.get("Cc")+formula.get("Cc"));
    				classElements.put("C", classElements.get("C")-formula.get("Cc"));
    			}
    		}
    	}
    	if (!labels.isEmpty()) 
    	{
    		Collections.sort(labels);
      	StringBuilder builder = new StringBuilder();
      	for (String label : labels) 
    		{
    			builder.append(label);
    		}
      	String labelID = builder.toString();
      	if (!labelElements.keySet().contains(labelID))
      	{
      		labelElements.put(labelID, classElements);
      	}
    	}
    }
  	return labelElements;
  }
	
	private Hashtable<String,Hashtable<String,Integer>> fillUnlabeledElements(Hashtable<String,Hashtable<String,Integer>> labelElements, 
  		int fattyAcidC, int dbs, String lipidClass)
  {
  	String label = "";
  	Hashtable<String,Integer> elements = new Hashtable<String,Integer>();
  	switch (lipidClass)
  	{
  		case "PC":
  			elements.put("C", fattyAcidC+8);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs);
  			elements.put("D", 0);
  			elements.put("O", 8);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "O-PC":
  			elements.put("C", fattyAcidC+8);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs+2);
  			elements.put("D", 0);
  			elements.put("O", 7);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "P-PC":
  			elements.put("C", fattyAcidC+8);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs);
  			elements.put("D", 0);
  			elements.put("O", 7);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "LPC":
  			elements.put("C", fattyAcidC+8);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs+2);
  			elements.put("D", 0);
  			elements.put("O", 7);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "PE":
  			elements.put("C", fattyAcidC+5);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs);
  			elements.put("D", 0);
  			elements.put("O", 8);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "O-PE":
  			elements.put("C", fattyAcidC+5);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs+2);
  			elements.put("D", 0);
  			elements.put("O", 7);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "P-PE":
  			elements.put("C", fattyAcidC+5);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs);
  			elements.put("D", 0);
  			elements.put("O", 7);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "LPE":
  			elements.put("C", fattyAcidC+5);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs+2);
  			elements.put("D", 0);
  			elements.put("O", 7);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "PS":
  			elements.put("C", fattyAcidC+6);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-2);
  			elements.put("D", 0);
  			elements.put("O", 10);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "LPS":
  			elements.put("C", fattyAcidC+6);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs);
  			elements.put("D", 0);
  			elements.put("O", 9);
  			elements.put("P", 1);
  			elements.put("N", 1);
  			break;
  		case "PI":
  			elements.put("C", fattyAcidC+9);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-3);
  			elements.put("D", 0);
  			elements.put("O", 13);
  			elements.put("P", 1);
  			elements.put("N", 0);
  			break;
  		case "LPI":
  			elements.put("C", fattyAcidC+9);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-1);
  			elements.put("D", 0);
  			elements.put("O", 12);
  			elements.put("P", 1);
  			elements.put("N", 0);
  			break;
  		case "PG":
  			elements.put("C", fattyAcidC+6);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-1);
  			elements.put("D", 0);
  			elements.put("O", 10);
  			elements.put("P", 1);
  			elements.put("N", 0);
  			break;
  		case "LPG":
  			elements.put("C", fattyAcidC+6);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs+1);
  			elements.put("D", 0);
  			elements.put("O", 9);
  			elements.put("P", 1);
  			elements.put("N", 0);
  			break;
  		case "SM":
  			elements.put("C", fattyAcidC+5);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs+3);
  			elements.put("D", 0);
  			elements.put("O", 6);
  			elements.put("P", 1);
  			elements.put("N", 2);
  			break;
  		case "Cer":
  			elements.put("C", fattyAcidC);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs+1);
  			elements.put("D", 0);
  			elements.put("O", 3);
  			elements.put("P", 0);
  			elements.put("N", 1);
  			break;
  		case "GM1": //done
  			elements.put("C", fattyAcidC+37);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-13);
  			elements.put("D", 0);
  			elements.put("O", 31);
  			elements.put("P", 0);
  			elements.put("N", 3);
  			break;
  		case "GM1-Ac": //done
  			elements.put("C", fattyAcidC+37+2);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-13-2);
  			elements.put("D", 0);
  			elements.put("O", 31+1);
  			elements.put("P", 0);
  			elements.put("N", 3);
  			break;
  		case "GM1-Fuc": //done
  			elements.put("C", fattyAcidC+37+12);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-13-4);
  			elements.put("D", 0);
  			elements.put("O", 31+9);
  			elements.put("P", 0);
  			elements.put("N", 3);
  			break;
  		case "GM2": //done
  			elements.put("C", fattyAcidC+31);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-11);
  			elements.put("D", 0);
  			elements.put("O", 26);
  			elements.put("P", 0);
  			elements.put("N", 3);
  			break;
  		case "GM3": //done
  			elements.put("C", fattyAcidC+23);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-8);
  			elements.put("D", 0);
  			elements.put("O", 21);
  			elements.put("P", 0);
  			elements.put("N", 2);
  			break;
  		case "GM4": //done
  			elements.put("C", fattyAcidC+17);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-8);
  			elements.put("D", 0);
  			elements.put("O", 16);
  			elements.put("P", 0);
  			elements.put("N", 2);
  			break;
  		case "GD1": //done
  			elements.put("C", fattyAcidC+48);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-18);
  			elements.put("D", 0);
  			elements.put("O", 39);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GD1-Ac": //done
  			elements.put("C", fattyAcidC+48+2);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-18-2);
  			elements.put("D", 0);
  			elements.put("O", 39+1);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GD1-Fuc": //done
  			elements.put("C", fattyAcidC+48+12);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-18-4);
  			elements.put("D", 0);
  			elements.put("O", 39+9);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GD1-Fuc-Ac": //done
  			elements.put("C", fattyAcidC+48+8);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-18-4);
  			elements.put("D", 0);
  			elements.put("O", 39+5);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GD2": //done
  			elements.put("C", fattyAcidC+42);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-16);
  			elements.put("D", 0);
  			elements.put("O", 34);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GD3": //done
  			elements.put("C", fattyAcidC+34);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-13);
  			elements.put("D", 0);
  			elements.put("O", 29);
  			elements.put("P", 0);
  			elements.put("N", 3);
  			break;
  		case "GD3-Ac": //done
  			elements.put("C", fattyAcidC+34+2);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-13-2);
  			elements.put("D", 0);
  			elements.put("O", 29+1);
  			elements.put("P", 0);
  			elements.put("N", 3);
  			break;
  		case "GT1": //done
  			elements.put("C", fattyAcidC+59);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-23);
  			elements.put("D", 0);
  			elements.put("O", 47);
  			elements.put("P", 0);
  			elements.put("N", 5);
  			break;
  		case "GT1-Ac": //done
  			elements.put("C", fattyAcidC+59+2);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-23-2);
  			elements.put("D", 0);
  			elements.put("O", 47+1);
  			elements.put("P", 0);
  			elements.put("N", 5);
  			break;
  		case "GT1-Fuc": //done
  			elements.put("C", fattyAcidC+59+12);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-23-4);
  			elements.put("D", 0);
  			elements.put("O", 47+9);
  			elements.put("P", 0);
  			elements.put("N", 5);
  			break;
  		case "GT2": //done
  			elements.put("C", fattyAcidC+53);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-21);
  			elements.put("D", 0);
  			elements.put("O", 42);
  			elements.put("P", 0);
  			elements.put("N", 5);
  			break;
  		case "GT3": //done
  			elements.put("C", fattyAcidC+45);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-18);
  			elements.put("D", 0);
  			elements.put("O", 37);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GT3-Ac": //done
  			elements.put("C", fattyAcidC+45+2);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-18-2);
  			elements.put("D", 0);
  			elements.put("O", 37+1);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GT3-Fuc": //done
  			elements.put("C", fattyAcidC+45+12);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-18-4);
  			elements.put("D", 0);
  			elements.put("O", 37+9);
  			elements.put("P", 0);
  			elements.put("N", 4);
  			break;
  		case "GQ1": //done
  			elements.put("C", fattyAcidC+70);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-28);
  			elements.put("D", 0);
  			elements.put("O", 55);
  			elements.put("P", 0);
  			elements.put("N", 6);
  			break;
  		case "GQ1-Ac": //done
  			elements.put("C", fattyAcidC+70+2);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-28-2);
  			elements.put("D", 0);
  			elements.put("O", 55+1);
  			elements.put("P", 0);
  			elements.put("N", 6);
  			break;
  		case "GQ1-Fuc": //done
  			elements.put("C", fattyAcidC+70+12);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-28-4);
  			elements.put("D", 0);
  			elements.put("O", 55+9);
  			elements.put("P", 0);
  			elements.put("N", 6);
  			break;
  		case "GP1": //done
  			elements.put("C", fattyAcidC+81);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-33);
  			elements.put("D", 0);
  			elements.put("O", 63);
  			elements.put("P", 0);
  			elements.put("N", 7);
  			break;
  		case "GH1": //done
  			elements.put("C", fattyAcidC+92);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-38);
  			elements.put("D", 0);
  			elements.put("O", 71);
  			elements.put("P", 0);
  			elements.put("N", 8);
  			break;
  		case "GS1": //done
  			elements.put("C", fattyAcidC+103);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-43);
  			elements.put("D", 0);
  			elements.put("O", 79);
  			elements.put("P", 0);
  			elements.put("N", 9);
  			break;
  		case "GS1-Fuc": //done
  			elements.put("C", fattyAcidC+103+12);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-43-4);
  			elements.put("D", 0);
  			elements.put("O", 79+9);
  			elements.put("P", 0);
  			elements.put("N", 9);
  			break;
  		case "GO1": //done
  			elements.put("C", fattyAcidC+114);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-48);
  			elements.put("D", 0);
  			elements.put("O", 87);
  			elements.put("P", 0);
  			elements.put("N", 10);
  			break;
  		case "LacCer": //done
  			elements.put("C", fattyAcidC+12);
  			elements.put("Cc", 0);
  			elements.put("H", elements.get("C")*2-2*dbs-3);
  			elements.put("D", 0);
  			elements.put("O", 13);
  			elements.put("P", 0);
  			elements.put("N", 1);
  			break;
  		default:
  			break;
  	}
  	labelElements.put(label, elements);
  	return labelElements;
  }
	
	private Vector<Vector<FattyAcidVO>> getPossCombis(int chains, int cAtoms, int dbs, Vector<FattyAcidVO> fas, Vector<FattyAcidVO> added){
    Vector<Vector<FattyAcidVO>> combis = new Vector<Vector<FattyAcidVO>>();
    for (FattyAcidVO fa : fas) {
      Vector<FattyAcidVO> toAdd = new Vector<FattyAcidVO>(added);
      toAdd.add(fa);
      if (chains>1) {
        Vector<Vector<FattyAcidVO>> combis2 = getPossCombis(chains-1, cAtoms, dbs, fas, toAdd);
        if (combis2.size()>0) {
          combis.addAll(combis2);
        }
      }else{
        if (checkTotalCDbsValid(toAdd,cAtoms,dbs)) {
          combis.add(toAdd);
        }
      }
    }
    return combis;
  }
	
	private boolean checkTotalCDbsValid(Vector<FattyAcidVO> toAdd, int cAtoms, int dbs) {
    int totC = 0;
    int totD = 0;
    for (FattyAcidVO fa : toAdd) {
      totC += fa.getcAtoms();
      totD += fa.getDoubleBonds();
    }
    return (totC==cAtoms && totD==dbs);
  }
}
