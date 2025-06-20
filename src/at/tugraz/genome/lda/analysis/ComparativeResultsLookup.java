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

import java.util.Vector;

import at.tugraz.genome.lda.vos.ResultAreaVO;

/**
 * the lookup to the comparative values over several LDA results (used in Statistical Analysis section)
 * 
 * @author Juergen Hartler
 *
 */
public interface ComparativeResultsLookup
{
  /**
   * returns the corresponding ResultAreaVO of a heat map
   * @param molGroup the analyte group
   * @param molName the analyte name (including the retention time if appropriate)
   * @param expName the name of the experiment
   * @return
   */
  public ResultAreaVO getResultAreaVO (String molGroup, String molName, String expName);
  
  /**
   * 
   * @return the experiment names in the sequence they are selected by the user
   */
  public Vector<String> getExpNamesInSequence();
}
