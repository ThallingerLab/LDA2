/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger 
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

package at.tugraz.genome.parsers;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.Vector;

import at.tugraz.genome.exception.LipidBLASTException;
import at.tugraz.genome.vos.LipidBLASTDetectionVO;
import at.tugraz.genome.vos.LipidBLASTIdentificationVO;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidBLASTParser
{
  
  private String fileName_;
  private final static int NAME_MS1_COLON = 0;
  private final static int NAME_ADDUCT_COLON = 1;
  private final static int NAME_MS2_COLON = 2;
  
  private Hashtable<String,Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>> results_;
  
  public LipidBLASTParser(String fileName){
    this.fileName_ = fileName;
  }
  
  public void parse() throws LipidBLASTException{
    results_  = new Hashtable<String,Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>>();
    File file = new File(fileName_);
    if (!file.exists()) throw new LipidBLASTException("The file \""+fileName_+"\" does not exist!");
    if (file.isDirectory()) throw new LipidBLASTException("The file \""+fileName_+"\" is a directory!");
    LineNumberReader reader = null;
    Vector<LipidBLASTDetectionVO> detectVOs = new Vector<LipidBLASTDetectionVO>();
    try{
      reader = new LineNumberReader(new FileReader(fileName_));
      String line;
      boolean headerFound = false;
      int scanNumColumn = -1;
      int fileColumn = -1;
      int precursorColumn = -1;
      int rankColumn = -1;
      int libraryColumn = -1;
      int libraryIdColumn = -1;
      int massColumn = -1;
      int deltaColumn = -1;
      int libMzColumn = -1;
      int scoreColumn = -1;
      int dotProductColumn = -1;
      int probabilityColumn = -1;
      int revDotColumn = -1;
      int peptideColumn = -1;
      int minimumAmountOfColumns = -1;
      int lineNumber = 0;
      while ((line=reader.readLine())!=null){
        lineNumber++;
        if (!headerFound){
          if (line.indexOf("Peptide")!=-1 && line.indexOf("Unknown")!=-1 && line.indexOf("Unknown")!=-1 && 
              line.indexOf("Rank")!=-1 && line.indexOf("Prob(%)")!=-1){
            headerFound = true;
            StringTokenizer tokenizer = new StringTokenizer(line,"\t");
            int columnCount = 0;
            while (tokenizer.hasMoreTokens()){
              String header = tokenizer.nextToken().trim();
              if (header.equalsIgnoreCase("Num")) scanNumColumn = columnCount;
              else if (header.equalsIgnoreCase("Unknown")) fileColumn = columnCount;
              else if (header.equalsIgnoreCase("Precursor m/z")) precursorColumn = columnCount;
              else if (header.equalsIgnoreCase("Rank")) rankColumn = columnCount;
              else if (header.equalsIgnoreCase("Library")) libraryColumn = columnCount;
              else if (header.equalsIgnoreCase("Id")) libraryIdColumn = columnCount;
              else if (header.equalsIgnoreCase("Mass")) massColumn = columnCount;
              else if (header.equalsIgnoreCase("Delta(m/z)")) deltaColumn = columnCount;
              else if (header.equalsIgnoreCase("Lib Precursor m/z")) libMzColumn = columnCount;
              else if (header.equalsIgnoreCase("Score")) scoreColumn = columnCount;
              else if (header.equalsIgnoreCase("Dot Product")) dotProductColumn = columnCount;
              else if (header.equalsIgnoreCase("Prob(%)")) probabilityColumn = columnCount;
              else if (header.equalsIgnoreCase("Rev-Dot")) revDotColumn = columnCount;
              else if (header.equalsIgnoreCase("Peptide")) peptideColumn = columnCount;
              columnCount++;
            }
            if (peptideColumn>minimumAmountOfColumns) minimumAmountOfColumns = peptideColumn;
            if (revDotColumn>minimumAmountOfColumns) minimumAmountOfColumns = revDotColumn;
            if (probabilityColumn>minimumAmountOfColumns) minimumAmountOfColumns = probabilityColumn;
            if (dotProductColumn>minimumAmountOfColumns) minimumAmountOfColumns = dotProductColumn;
            if (scoreColumn>minimumAmountOfColumns) minimumAmountOfColumns = scoreColumn;
            if (libMzColumn>minimumAmountOfColumns) minimumAmountOfColumns = libMzColumn;
            if (deltaColumn>minimumAmountOfColumns) minimumAmountOfColumns = deltaColumn;
            if (massColumn>minimumAmountOfColumns) minimumAmountOfColumns = massColumn;
            if (libraryIdColumn>minimumAmountOfColumns) minimumAmountOfColumns = libraryIdColumn;
            if (libraryColumn>minimumAmountOfColumns) minimumAmountOfColumns = libraryColumn;
            if (rankColumn>minimumAmountOfColumns) minimumAmountOfColumns = rankColumn;
            if (precursorColumn>minimumAmountOfColumns) minimumAmountOfColumns = precursorColumn;
            if (fileColumn>minimumAmountOfColumns) minimumAmountOfColumns = fileColumn;
            if (scanNumColumn>minimumAmountOfColumns) minimumAmountOfColumns = scanNumColumn;            
          }
        }else{
          StringTokenizer tokenizer = new StringTokenizer(line,"\t");
          if (tokenizer.countTokens()<minimumAmountOfColumns) continue;
          String ms1Name = null;
          String adduct = null;
          String ms2Name = null;
          String retentionTime = null;
          int scanNum = -1;
          String fileName = null;
          String precursorMass = null;
          int rank = -1;
          String library = null;
          long libraryId = -1;
          String mass = null;
          String deltaMz = null;
          String libMz = null;
          int score = -1;
          int dotProduct = -1;
          double probability = -1;
          int revDot = -1;
          int columnCount = 0;
          while (tokenizer.hasMoreTokens()){
            String item = tokenizer.nextToken();
            if (columnCount==peptideColumn){
              StringTokenizer tok2 = new StringTokenizer(item,";");
              if (tok2.countTokens()<3) System.out.println("Attention: the Peptide column \""+item+"\" does not contain enough information");
              int colonCount = 0;
              while (tok2.hasMoreTokens()){
                String token = tok2.nextToken().trim();
                if (colonCount==NAME_MS1_COLON) ms1Name = token;
                else if (colonCount==NAME_ADDUCT_COLON) adduct = token;
                else if (colonCount==NAME_MS2_COLON) ms2Name = token;
                colonCount++;
              }
            } else if (columnCount==fileColumn){
              if (item.indexOf(" RT:")!=-1){
                fileName = item.substring(0,item.indexOf(" RT:")).trim();
                retentionTime = item.substring(item.indexOf(" RT:")+" RT:".length()).trim();
              }
            } else if (columnCount==scanNumColumn){
              try{
                scanNum = Integer.parseInt(item);
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A scan number must be integer format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==precursorColumn){
              try{
                Double.parseDouble(item);
                precursorMass = item;
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Precursor m/z\" must be double format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==rankColumn){
              try{
                rank = Integer.parseInt(item);
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Rank\" must be integer format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==libraryColumn){
              library = item;
            } else if (columnCount==libraryIdColumn){
              try{
                libraryId = Long.parseLong(item);
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A library ID must be long format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==massColumn){
              try{
                Double.parseDouble(item);
                mass = item;
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Mass\" must be double format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==deltaColumn){
              try{
                Double.parseDouble(item);
                deltaMz = item;
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Delta(m/z)\" must be double format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==libMzColumn){
              try{
                Double.parseDouble(item);
                libMz = item;
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Lib Precursor m/z\" must be double format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==scoreColumn){
              try{
                score = Integer.parseInt(item);
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Score\" must be integer format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==dotProductColumn){
              try{
                dotProduct = Integer.parseInt(item);
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Dot Product\" must be integer format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==probabilityColumn){
              try{
                probability = Double.parseDouble(item);
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Prob(%)\" must be double format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            } else if (columnCount==revDotColumn){
              try{
                revDot = Integer.parseInt(item);
              } catch (NumberFormatException nfx){
                throw new LipidBLASTException("A \"Rev-Dot\" must be integer format! The one ("+item+") at line number "+lineNumber+" is not!");
              }
            }
            columnCount++;
          }
          if (ms1Name!=null && ms1Name.length()>0 && adduct!=null && adduct.length()>0 && ms2Name!=null && ms2Name.length()>0 &&
              retentionTime!=null && retentionTime.length()>0 && rank>-1 && probability>0){
            LipidBLASTDetectionVO detectVO = new LipidBLASTDetectionVO(ms1Name, adduct, ms2Name, retentionTime, scanNum,
                fileName, precursorMass, rank, library, libraryId, mass, deltaMz, libMz, score, dotProduct, probability, revDot);
            detectVOs.add(detectVO);
          }
        }
      }
    } catch (IOException iox){
      throw new LipidBLASTException(iox);
    } finally {
      if (reader!=null){
        try {
          reader.close();
        }
        catch (IOException e){}
      }  
    }
    //build the result hash
    for (LipidBLASTDetectionVO detectVO : detectVOs){
      Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>> classResults = new Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>();
      if (results_.containsKey(detectVO.getLipidClass())) classResults = results_.get(detectVO.getLipidClass());
      Hashtable<String,LipidBLASTIdentificationVO> analyteResults = new Hashtable<String,LipidBLASTIdentificationVO>();
      if (classResults.containsKey(detectVO.getMs1Name())) analyteResults = classResults.get(detectVO.getMs1Name());
      LipidBLASTIdentificationVO idVO = new LipidBLASTIdentificationVO();
      if (analyteResults.containsKey(detectVO.getAdduct())) idVO = analyteResults.get(detectVO.getAdduct());
      idVO.addDetection(detectVO);
      analyteResults.put(detectVO.getAdduct(), idVO);
      classResults.put(detectVO.getMs1Name(), analyteResults);
      results_.put(detectVO.getLipidClass(), classResults);
    }
  }

  public Hashtable<String,Hashtable<String,Hashtable<String,LipidBLASTIdentificationVO>>> getResults_()
  {
    return results_;
  }
  
  
}
