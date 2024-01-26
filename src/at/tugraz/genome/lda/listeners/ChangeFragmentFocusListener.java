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

import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import at.tugraz.genome.lda.swing.RuleDefinitionInterface;

/**
 * Listener if the focus is lost.
 * @author Andreas Ziegl
 *
 */
public class ChangeFragmentFocusListener implements FocusListener 
{
	/** The RDI self */
	private RuleDefinitionInterface rDI_;
	
	/** The position of the fragment in the table */
  private int position_; 
  
  /** The type of the fragment info (Formula, Charge,...) */
  private int type_; 
  
  /** Head or chain fragment */
  private int headOrChain_;

  /**
   * The constuctor of the listener with alls params for the checkFragmentWithErrors function
   * @param rDI
   * @param position
   * @param type
   * @param headOrChain
   */
  public ChangeFragmentFocusListener(RuleDefinitionInterface rDI, int position, int type, int headOrChain) 
  {
      this.rDI_ = rDI;
      this.position_ = position;
      this.type_ = type; 
      this.headOrChain_ = headOrChain;
  }

	@Override
	public void focusLost(FocusEvent e)
	{
		rDI_.checkFragmentWithErrors(position_, type_, headOrChain_);		
	}

	@Override
	public void focusGained(FocusEvent e)
	{		
	}

}
