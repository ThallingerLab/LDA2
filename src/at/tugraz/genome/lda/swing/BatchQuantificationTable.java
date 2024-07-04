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

package at.tugraz.genome.lda.swing;

import java.awt.Rectangle;

import javax.swing.JTable;
import javax.swing.JViewport;
import javax.swing.table.TableModel;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class BatchQuantificationTable extends JTable
{

  /**
   * 
   */
  private static final long serialVersionUID = 3369368275154118727L;

  public BatchQuantificationTable(TableModel model)
  {
    super(model);
  }

  /**
     * Scrolls the cell (rowIndex, colIndex) so that it is visible at the center
     * of viewport. Assumes table is contained in a JScrollPane.
     * 
     * @param rowIndex
     * @param colIndex
     * @author javaalmanac.com
     */
  public void scrollToCenter(int rowIndex, int colIndex)
  {
    if (!(getParent() instanceof JViewport)) {
      return;
    }
    JViewport viewport = (JViewport) getParent();

    // This rectangle is relative to the table where the
    // northwest corner of cell (0,0) is always (0,0).
    Rectangle rect = getCellRect(rowIndex, colIndex, true);

    // The location of the view relative to the table
    Rectangle viewRect = viewport.getViewRect();

    // Translate the cell location so that it is relative
    // to the view, assuming the northwest corner of the
    // view is (0,0).
    rect.setLocation(rect.x - viewRect.x, rect.y - viewRect.y);

    // Calculate location of rect if it were at the center of view
    int centerX = (viewRect.width - rect.width) / 2;
    int centerY = (viewRect.height - rect.height) / 2;

    // Fake the location of the cell so that scrollRectToVisible
    // will move the cell to the center
    if (rect.x < centerX) {
      centerX = -centerX;
    }
    if (rect.y < centerY) {
      centerY = -centerY;
    }
    rect.translate(centerX, centerY);

    // Scroll the area into view.
    try {
      viewport.scrollRectToVisible(rect);
    }
    catch (Exception ex) {
      ex.printStackTrace();
    }
  }

}
