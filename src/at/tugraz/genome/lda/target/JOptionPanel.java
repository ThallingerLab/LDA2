/* 
 * This file is part of Lipid Data Analyzer
 * Lipid Data Analyzer - Automated annotation of lipid species and their molecular structures in high-throughput data from tandem mass spectrometry
 * Copyright (c) 2023 Juergen Hartler, Andreas Ziegl, Gerhard G. Thallinger, Leonida M. Lamp
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

package at.tugraz.genome.lda.target;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JPanel;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;

/**
 * 
 * @author Leonida M. Lamp
 *
 */
public class JOptionPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  public final static String ICON_LCCL_PATH = "/images/Logo.png"; 
  
  protected JDefaultComponents defaultComponents_;
  private String panelTitle_;
  
  protected JOptionPanel(JDefaultComponents defaultComponents) 
  {
  	this.defaultComponents_ = defaultComponents;
    this.panelTitle_ = null;
  }
  
  protected JOptionPanel(JDefaultComponents defaultComponents, String title)
  {
  	this.defaultComponents_ = defaultComponents;
    this.panelTitle_ = title;
  }
  
  protected void update() 
  {
  }
  
  protected void next() 
  {
    goNext();
  }
  
  protected void back() 
  {
    goBack();
  }
  
  protected boolean goNext() 
  {
    if (defaultComponents_.getOptionPanelList().size() > defaultComponents_.getCurrentIndex()+1 ) 
    {
      defaultComponents_.setCurrentIndex(defaultComponents_.getCurrentIndex()+1);
      defaultComponents_.updateComponents();
      return true;
    } 
    else 
    {
      return false;
    }
  }
  
  protected boolean goBack() 
  {
    if (defaultComponents_.getCurrentIndex()-1 >= 0) 
    {
      defaultComponents_.setCurrentIndex(defaultComponents_.getCurrentIndex()-1);
      defaultComponents_.updateComponents();
      return true;
    } 
    else 
    {
      return false;
    }
  }
  
  protected JDefaultComponents getDefaultComponents()
  {
    return defaultComponents_;
  }
  
  protected void setDefaultComponents(JDefaultComponents aDefaultComponents)
  {
    defaultComponents_ = aDefaultComponents;
  }
  
  protected String getPanelTitle() 
  {
    return panelTitle_;
  }
  
  protected void setPanelTitle(String title) 
  {
    panelTitle_ = title;
  }
  
  protected ImageIcon getDefaultImageIcon(String path,
      String description) 
  {
		java.net.URL imgURL = JOptionPanel.class.getResource(path);
		if (imgURL != null) 
		{
			return new ImageIcon(imgURL, description);
		} 
		else 
		{
			System.err.println("Couldn't find file: " + path);
			return null;
		}
  }	
  
  protected void switchPanel(int panelIndex) 
  {
  	defaultComponents_.setCurrentIndex(panelIndex);
  	defaultComponents_.updateComponents();
  }
  
  public static TitledBorder getTitledPanelBorder(String title)
  {
  	Border raisedbevel = BorderFactory.createRaisedBevelBorder();
  	Border loweredbevel = BorderFactory.createLoweredBevelBorder();
  	Border compound = BorderFactory.createCompoundBorder(raisedbevel, loweredbevel);
  	return BorderFactory.createTitledBorder(compound, title);
  }
  
  
}
