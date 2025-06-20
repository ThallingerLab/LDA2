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

package at.tugraz.genome.lda.analysis;

import java.util.LinkedHashMap;
import java.util.Vector;

/**
 * 
 * @author Juergen Hartler
 *
 */
public interface SampleLookup
{
  public String getDisplayName(String sampleName);
  public void setDisplayName(String sampleName, String displayName);
  
  /**
   * 
   * @return a sorted hash map containing the full file paths of the LDA-results; key: abbreviated experiment name; value: full file path
   */
  public LinkedHashMap<String,String> getSampleResultFullPaths();
  
  /**
   * 
   * @return a sorted hash map containing the sample groups as key, and the abbreviated experiment names belonging to this group as values
   */
  public LinkedHashMap<String,Vector<String>> getSamplesOfGroups();
  
  /**
   * 
   * @return the lookup to the comparative values; i.e. the (RT-grouped) result area VOs, and the original experiment names
   */
  public ComparativeResultsLookup getComparativeResultsLookup();
  
}
