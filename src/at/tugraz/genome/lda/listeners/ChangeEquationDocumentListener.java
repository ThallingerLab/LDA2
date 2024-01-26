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

package at.tugraz.genome.lda.listeners;

import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import at.tugraz.genome.lda.swing.RuleDefinitionInterface;

/**
 * Listener check the added equations while the user is entering
 * Also paints a new spectra
 * @author Andreas Ziegl
 *
 */
public class ChangeEquationDocumentListener implements DocumentListener 
{
	/** The RDI self */
  private RuleDefinitionInterface rDI_;
  
  /** The position in the table */
  private int position_;
  
  /** Head or chain fragment */
  private int headOrChainOrPostion_;

    /**
     * The constructor of the listener.
     * @param rDI
     * @param position
     * @param headOrChainPosition
     */
    public ChangeEquationDocumentListener(RuleDefinitionInterface rDI, int position, int headOrChainOrPostion) 
    {
        this.rDI_ = rDI;
        this.position_ = position;
        this.headOrChainOrPostion_ = headOrChainOrPostion;
    }

    @Override
		public void insertUpdate(DocumentEvent e)
    {
   	  rDI_.refreshEquationFieldToVerifie(position_, headOrChainOrPostion_);
    }
    
    @Override
    public void removeUpdate(DocumentEvent e)
    {
   	  rDI_.refreshEquationFieldToVerifie(position_, headOrChainOrPostion_);
    }  
    
    @Override
    public void changedUpdate(DocumentEvent e)
    {
   	  rDI_.refreshEquationFieldToVerifie(position_, headOrChainOrPostion_);
    }  
}

