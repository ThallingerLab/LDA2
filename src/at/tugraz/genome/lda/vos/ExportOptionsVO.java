/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2017 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp 
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

package at.tugraz.genome.lda.vos;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class ExportOptionsVO
{
  
  private int exportType_;
  private String sdValue_;
  private boolean analyteInColumn_;
  private boolean exportRT_;
  private boolean exportRTDev_;
  private boolean exportDoubleBondPositions_;
  private int commaPositions_;
  /** which type of species shall be exported*/
  private short speciesType_;
  
  public static final int EXPORT_NO_DEVIATION = 0;
  public static final int EXPORT_SD_DEVIATION = 1;
  public static final int EXPORT_SD_ERROR = 2;
  public static final int EXPORT_SD_DEV_AND_ERROR = 3;
  
  
  /**
   * Constructor for object containing the user-selected export options
   * @param exportType which values shall be exported (see ExportOptionsVO.EXPORT_ ...)
   * @param sdValue the multiplicator for the standard deviation
   * @param analyteInColumn true when analytes are in the columns, false when they are in the rows
   * @param exportRT true when retention time values shall be exported
   * @param exportRTDev true when standard deviation values of retention times shall be exported
   * @param exportDoubleBondPositions true when double bond positions shall be exported
   * @param commaPositions the amounts of decimal places after the comma
   * @param speciesType structural level of data export (lipid species, chain level, position level - for details see LipidomicsConstants.EXPORT_ANALYTE_TYPE)
   */
  public ExportOptionsVO(int exportType,String sdValue, boolean analyteInColumn, boolean exportRT, boolean exportRTDev, boolean exportDoubleBondPositions, 
      int commaPositions, short speciesType){
    exportType_ = exportType;
    sdValue_ = sdValue;
    analyteInColumn_ = analyteInColumn;
    exportRT_ = exportRT;
    exportRTDev_ = exportRTDev;
    exportDoubleBondPositions_ = exportDoubleBondPositions;
    commaPositions_ = commaPositions;
    speciesType_ = speciesType;
  }



  public int getExportType()
  {
    return exportType_;
  }

  public String getSdValue()
  {
    return sdValue_;
  }
  

  public boolean isAnalyteInColumn()
  {
    return analyteInColumn_;
  }


  public boolean isExportRT()
  {
    return exportRT_;
  }


  public boolean isExportRTDev()
  {
    return exportRTDev_;
  }
  
  public boolean isExportDoubleBondPositions()
  {
    return exportDoubleBondPositions_;
  }



  public int getCommaPositions()
  {
    return commaPositions_;
  }


  /**
   * 
   * @return which type of species shall be exported
   */
  public short getSpeciesType()
  {
    return speciesType_;
  }

  
  
}
