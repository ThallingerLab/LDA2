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
  private int commaPositions_;
  
  public static final int EXPORT_NO_DEVIATION = 0;
  public static final int EXPORT_SD_DEVIATION = 1;
  public static final int EXPORT_SD_ERROR = 2;
  public static final int EXPORT_SD_DEV_AND_ERROR = 3;
  
  public ExportOptionsVO(int exportType,String sdValue, boolean analyteInColumn, boolean exportRT, boolean exportRTDev, int commaPositions){
    exportType_ = exportType;
    sdValue_ = sdValue;
    analyteInColumn_ = analyteInColumn;
    exportRT_ = exportRT;
    exportRTDev_ = exportRTDev;
    commaPositions_ = commaPositions;
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



  public int getCommaPositions()
  {
    return commaPositions_;
  }
  
  
  
}
