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

package at.tugraz.genome.lda.swing;

import java.awt.Color;
import java.awt.Component;

import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.table.DefaultTableCellRenderer;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsTableCellRenderer extends DefaultTableCellRenderer
{
  private static final long serialVersionUID = 393014627701796523L;
  
  public final static Color BRIGHT_GREEN = new Color(200,255,200);
  


  public Component getTableCellRendererComponent(JTable table, Object
      value, boolean isSelected,
      boolean hasFocus, int
      row, int column) {
      JLabel renderedLabel = (JLabel) super.getTableCellRendererComponent(
      table, value, isSelected, hasFocus, row, column);
      // the name is aligned left and the area right 
      if (column == LipidomicsTableModel.COLUMN_AREA) renderedLabel.setHorizontalAlignment(JLabel.RIGHT);
      else renderedLabel.setHorizontalAlignment(JLabel.LEFT);
      LipidomicsTableModel model = (LipidomicsTableModel)table.getModel();
      if (!isSelected){
        
        if (model.isPercentalSplitInstance(row)) renderedLabel.setBackground(Color.ORANGE);
        else if (model.isSplitInstance(row)) renderedLabel.setBackground(Color.YELLOW);
        // data containing MSn evidence is colored in green
        else if (model.hasMS2Evidence(row)) renderedLabel.setBackground(BRIGHT_GREEN);
        else renderedLabel.setBackground(Color.WHITE);
        if (model.hasOmegaEvidenceConflict(row)) renderedLabel.setForeground(Color.RED);
        else renderedLabel.setForeground(Color.BLACK);
        
      }
      return renderedLabel;
  }
  

}
