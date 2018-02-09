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

import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;
import javax.swing.SwingUtilities;
import javax.swing.table.TableCellRenderer;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.analysis.AnalyteAddRemoveListener;
import at.tugraz.genome.lda.quantification.LipidParameterSet;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class LipidomicsJTable extends JTable implements ActionListener
{
  private static final long serialVersionUID = 1924361162765772777L;
  
  public final static int ORDER_TYPE_AS_IS = 0;
  public final static int ORDER_TYPE_MZ = 1;
  public final static int ORDER_TYPE_INTENSITY = 2;
  
  private TableCellRenderer renderer;
  
  private JPopupMenu addItemPopup_;
  private AnalyteAddRemoveListener parentListener_;
  
  private boolean showRt_;
  
  private boolean controlDown_;
  private boolean shiftDown_;
  private int lastSelectedIndex_;
  
  private final static String LABEL_EDIT_RULES = "Edit MSn rule";
  /** display text for recalculating an MSn identification of a hit*/
  private final static String LABEL_RECALCULATE_MSN = "Recalculate MSn";
  /** display text for editing the retention time of a hit*/
  private final static String LABEL_EDIT_RT = "Edit Rt";

  public LipidomicsJTable(LipidomicsTableModel model,TableCellRenderer renderer,boolean showMs2,int orderType,boolean showRtInAddDialog,
      AnalyteAddRemoveListener parentListener)  {
//    super(rowData,columnNames);
    super(model);
    parentListener_ = parentListener;
    this.renderer = renderer;
    this.showRt_ = showRtInAddDialog;
    addItemPopup_ = new JPopupMenu("Apply Peak-RT to double peaks");
    JMenuItem item = new JMenuItem("Add analyte before");
    item.addActionListener(this);
    item.setVisible(orderType==ORDER_TYPE_AS_IS);
    addItemPopup_.add(item);
    item = new JMenuItem("Add analyte after");
    item.addActionListener(this);
    item.setVisible(orderType==ORDER_TYPE_AS_IS);
    addItemPopup_.add(item);
    item = new JMenuItem("Delete analyte");
    item.addActionListener(this);
    addItemPopup_.add(item);
    item = new JMenuItem("Show MS/MS");
    item.addActionListener(this);
    item.setVisible(showMs2);
    addItemPopup_.add(item);

    item = new JMenuItem("Sort by order");
    item.addActionListener(this);
    item.setVisible(orderType!=ORDER_TYPE_AS_IS);
    addItemPopup_.add(item);
    item = new JMenuItem("Sort by mass");
    item.addActionListener(this);
    item.setVisible(orderType!=ORDER_TYPE_MZ);
    addItemPopup_.add(item);
    item = new JMenuItem("Sort by intensity");
    item.addActionListener(this);
    item.setVisible(orderType!=ORDER_TYPE_INTENSITY);
    addItemPopup_.add(item);
    
    //Edit MSn rule
    item = new JMenuItem(LABEL_EDIT_RULES);
    item.addActionListener(this);
    addItemPopup_.add(item);

    //Recalculate an MSn identification
    if (!model.isShowMSn()){
      item = new JMenuItem(LABEL_RECALCULATE_MSN);
      item.addActionListener(this);
      addItemPopup_.add(item);
      
      item = new JMenuItem(LABEL_EDIT_RT);
      item.addActionListener(this);
      addItemPopup_.add(item);
    }
    
    this.add(addItemPopup_);
    controlDown_ = false;
    shiftDown_ = false;
    lastSelectedIndex_ = 0;
    
    this.addMouseListener( new MouseAdapter(){
      public void mouseClicked( MouseEvent e ){
        // Left mouse click
        if ( SwingUtilities.isLeftMouseButton( e ) || SwingUtilities.isRightMouseButton( e ) ){
          Point p = e.getPoint();
          // get the row index that contains that coordinate
          ListSelectionModel model = getSelectionModel();
          int rowNumber = rowAtPoint( p );
          if (e.isControlDown() || e.isShiftDown()){
            if (e.isControlDown()) controlDown_ = true;
            if (e.isShiftDown()) shiftDown_ = true;
            if (SwingUtilities.isRightMouseButton( e )){
              if (e.isShiftDown()){
                model.setAnchorSelectionIndex(lastSelectedIndex_);
                model.setLeadSelectionIndex(rowNumber);
              }
              if (e.isControlDown()) model.addSelectionInterval(rowNumber,rowNumber);
              addItemPopup_.show(e.getComponent(), e.getX(), e.getY());
            }
            lastSelectedIndex_ = rowNumber;
          }else{
            controlDown_ = false;
            shiftDown_ = false;
            parentListener_.listSelectionChanged(rowNumber);
            if (SwingUtilities.isRightMouseButton( e )){
              model.setSelectionInterval( rowNumber, rowNumber );
              addItemPopup_.show(e.getComponent(), e.getX(), e.getY());              
            }
          }
        }
      }
    });
    this.addKeyListener( new KeyAdapter(){
      public void keyReleased(KeyEvent e) {
        if ((controlDown_&& !e.isControlDown()) || (shiftDown_&& !e.isShiftDown())){
          controlDown_ = false;
          shiftDown_ = false;
          parentListener_.listSelectionChanged(lastSelectedIndex_);
        }
      }  
    });
    this.setToolTipText(TooltipTexts.DISPLAY_SELECTION_TABLE);
  }

  public TableCellRenderer getCellRenderer(int row, int column) {
    return renderer;
  }

  public void actionPerformed(ActionEvent e)
  {
    this.actionPerformed(e.getActionCommand());
  }

  private void actionPerformed(String actionCommand) 
  {
    if (actionCommand.equalsIgnoreCase("Add analyte before")){
      int positionToAddNewAnalyte = getSelectionModel().getLeadSelectionIndex();
      LipidParameterSet params = parentListener_.getAnalyteInTableAtPosition(positionToAddNewAnalyte);
      String displayString = params.getNameString();
      if (params.getModificationName()!=null&&params.getModificationName().length()>0)displayString+="_"+params.getModificationName();
      new AddAnalyteDialog(new JFrame(),"Enter new analyte", "Add a new analyte before "+displayString,params.getNameString(),params.Mz[0],params.getAnalyteFormula(),params.getModificationName(),params.getModificationFormula(),positionToAddNewAnalyte,showRt_,parentListener_);
    } else if (actionCommand.equalsIgnoreCase("Add analyte after")){
      int positionToAddNewAnalyte = getSelectionModel().getLeadSelectionIndex()+1;
      LipidParameterSet params = parentListener_.getAnalyteInTableAtPosition(positionToAddNewAnalyte-1);
      String displayString = params.getNameString();
      if (params.getModificationName()!=null&&params.getModificationName().length()>0)displayString+="_"+params.getModificationName();
      new AddAnalyteDialog(new JFrame(),"Enter new analyte", "Add a new analyte after "+displayString,params.getNameString(),params.Mz[0],params.getAnalyteFormula(),params.getModificationName(),params.getModificationFormula(),positionToAddNewAnalyte,showRt_,parentListener_);
    } else if (actionCommand.equalsIgnoreCase("Delete analyte")){
      int[] indices = getSelectedRows();
      String displayString = "\n";
      int count = 0;
      for (int ind : indices){
        count++;
        LipidParameterSet params = parentListener_.getAnalyteInTableAtPosition(ind);
        displayString += params.getNameString();
        if (params.getModificationName()!=null&&params.getModificationName().length()>0)displayString+="_"+params.getModificationName();
        displayString += ", ";
        if (count%8 == 0) displayString += "\n";
      }
      if (displayString.endsWith("\n")) displayString = displayString.substring(0,displayString.length()-1);
      if (displayString.length()>1) displayString = displayString.substring(0,displayString.length()-2);
      if (JOptionPane.showConfirmDialog(this, "Do you really want to delete "+displayString+"?") == JOptionPane.YES_OPTION){
        parentListener_.removeAnalyte(indices);
      }
    } else if (actionCommand.equalsIgnoreCase("Show MS/MS")){
      int position = getSelectionModel().getLeadSelectionIndex();
      parentListener_.showMs2(position);
    } else if (actionCommand.equalsIgnoreCase("Sort by order")){
      parentListener_.changeListSorting(ORDER_TYPE_AS_IS);      
    } else if (actionCommand.equalsIgnoreCase("Sort by mass")){
      parentListener_.changeListSorting(ORDER_TYPE_MZ);
    } else if (actionCommand.equalsIgnoreCase("Sort by intensity")){
      parentListener_.changeListSorting(ORDER_TYPE_INTENSITY);
    } else if (actionCommand.equalsIgnoreCase(LABEL_EDIT_RULES)){
      int position = getSelectionModel().getLeadSelectionIndex();
      parentListener_.newRule(position);     
    } else if (actionCommand.equalsIgnoreCase(LABEL_RECALCULATE_MSN)){
      int position = getSelectionModel().getLeadSelectionIndex();
      parentListener_.recalculateMSn(position);     
    } else if (actionCommand.equalsIgnoreCase(LABEL_EDIT_RT)){
      int position = getSelectionModel().getLeadSelectionIndex();
      parentListener_.editRt(position);     
    }
  }
  
  public String getSumLipidNameAt(int rowIndex){
    return ((LipidomicsTableModel)super.getModel()).getSumLipidNameAt(rowIndex);
  }
  
}
