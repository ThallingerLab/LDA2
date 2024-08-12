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

package at.tugraz.genome.lda.verifier;

import javax.swing.InputVerifier;
import javax.swing.JComponent;

import at.tugraz.genome.lda.swing.rdi.GeneralSettingsPanel;

/**
 * For Verifying any of the input components of the "General" tab (GSP = general settings panel) of the rule definition interface (RDI)
 * @author Juergen Hartler
 *
 */
public class GeneralSettingsVerifier extends InputVerifier
{

  /** the GSP containing the method for the check*/
  private GeneralSettingsPanel gsp_;
  
  /**
   * constructor taking a GSP for performing the check
   * @param rdi the GSP containing the method for the check
   */
  public GeneralSettingsVerifier(GeneralSettingsPanel rdi){
    gsp_ = rdi;
  }
  
  public boolean verify(JComponent input)
  {
    // if the updateGeneralEntries() is successful it returns true; otherwise the input is not OK
    return gsp_.updateGeneralEntries();
  }

}
