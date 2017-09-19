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

import java.io.File;

import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class RawQuantificationPairVO
{
  private File rawFile;
  private File quantFile;
  /** is the raw file a ABSciex wiff file*/
  private boolean fromWiff_;
  private String status;
  
  public File getRawFile()
  {
    return rawFile;
  }

  public File getQuantFile()
  {
    return quantFile;
  }

  public RawQuantificationPairVO(File rawFile, File quantFile)
  {
    super();
    this.rawFile = rawFile;
    this.quantFile = quantFile;
    this.fromWiff_ = false;
  }
  
  public RawQuantificationPairVO(File rawFile, File quantFile, boolean isFromWiff)
  {
    this(rawFile,quantFile);
    this.fromWiff_ = isFromWiff;
  }
  
  public String getRawFileName(){
    return StaticUtils.extractFileName(rawFile.getAbsolutePath());
  }
  
  public String getQuantFileName(){
    return StaticUtils.extractFileName(quantFile.getAbsolutePath());
  }

  public String getStatus()
  {
    return status;
  }

  public void setStatus(String status)
  {
    this.status = status;
  }

  /**
   * 
   * @return if the raw file is a ABSciex wiff file
   */
  public boolean isFromWiff()
  {
    return fromWiff_;
  }
  
}
