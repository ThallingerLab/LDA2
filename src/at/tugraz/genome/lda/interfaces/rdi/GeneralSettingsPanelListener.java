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

package at.tugraz.genome.lda.interfaces.rdi;

/**
 * Callback listener for the panel for entering the general settings
 * @author Juergen Hartler
 *
 */
public interface GeneralSettingsPanelListener
{
  /** stores the general settings when a new class was entered - required for entering the fragmentation rules*/
  public boolean performStorageOfInitialGeneralSettings();
  /** shows the other tabs that depend on the general settings*/
  public void refreshGeneralSettingsDependantDisplay();
  /** asks if the tabs which depend on the general settings are currently visible*/
  public boolean areGeneralSettingsDependingFieldsVisible();
  /** returns the name of the lipid that is used for defining the rules*/
  public String getAnalyteName();
  /** returns the lipid class name for which the rules are currently defined*/
  public String getLipidClassName();
  /** refreshes all information that depends on the general settings*/
  public void updateGeneralSettings();
}
