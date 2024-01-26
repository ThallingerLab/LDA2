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

import java.util.ArrayList;

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

import at.tugraz.genome.lda.Settings;
import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.analysis.AnalyteAddRemoveListener;
import at.tugraz.genome.lda.msn.LipidomicsMSnSet;
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
  
  /** a popup menu item for activating the MS2 view*/
  private JMenuItem showMS2_;
  /** a popup menu item for activating the MS1 view*/
  private JMenuItem showMS1_;
  
  private boolean showRt_;
  
  private boolean controlDown_;
  private boolean shiftDown_;
  private int lastSelectedIndex_;

  private final static String LABEL_ADD_ANALYTE_BEFORE = "Add analyte before";
  
  private final static String LABEL_ADD_ANALYTE_AFTER = "Add analyte after";
  
  private final static String LABEL_DELETE_ANALYTE = "Delete analyte";
  
  private final static String LABEL_SHOW_MSMS = "Show MS/MS";
  
  private final static String LABEL_SHOW_MS1 = "Show MS1";
  
  private final static String LABEL_SORT_BY_ORDER = "Sort by order";
  
  private final static String LABEL_SORT_BY_MASS = "Sort by mass";
  
  private final static String LABEL_SORT_BY_INTENSITY = "Sort by intensity";
  
  private final static String LABEL_EDIT_RULES = "Edit MSn rule";
  /** display text for recalculating an MSn identification of a hit*/
  private final static String LABEL_RECALCULATE_MSN = "Recalculate MSn";
  /** display text for editing the retention time of a hit*/
  private final static String LABEL_EDIT_RT = "Edit Rt";
  
  private final static String LABEL_EDIT_OMEGA_ASSIGNMENT = "Edit \u03C9-DB Assignment";
  
  private final static String LABEL_DELETE_MOLECULAR_SPECIES = "Delete molecular species";
  /** do these LipidParameterSets contain any OH information*/
  private boolean hasOh_;

  public LipidomicsJTable(LipidomicsTableModel model,TableCellRenderer renderer,boolean showMs2,int orderType,boolean showRtInAddDialog,
      AnalyteAddRemoveListener parentListener)  {
//    super(rowData,columnNames);
    super(model);
    parentListener_ = parentListener;
    this.renderer = renderer;
    this.showRt_ = showRtInAddDialog;
    addItemPopup_ = new JPopupMenu("Apply Peak-RT to double peaks");
    
    JMenuItem item = new JMenuItem(LABEL_ADD_ANALYTE_BEFORE);
    item.addActionListener(this);
    item.setVisible(orderType==ORDER_TYPE_AS_IS);
    addItemPopup_.add(item);
    
    item = new JMenuItem(LABEL_ADD_ANALYTE_AFTER);
    item.addActionListener(this);
    item.setVisible(orderType==ORDER_TYPE_AS_IS);
    addItemPopup_.add(item);
    
    item = new JMenuItem(LABEL_DELETE_ANALYTE);
    item.addActionListener(this);
    addItemPopup_.add(item);
    
    if (model.isShowMSn()) {
      item = new JMenuItem(LABEL_DELETE_MOLECULAR_SPECIES);
      item.addActionListener(this);
      addItemPopup_.add(item);
    }
    
    showMS2_ = new JMenuItem(LABEL_SHOW_MSMS);
    showMS2_.addActionListener(this);
    showMS2_.setVisible(showMs2);
    addItemPopup_.add(showMS2_);
    
    showMS1_ = new JMenuItem(LABEL_SHOW_MS1);
    showMS1_.addActionListener(this);
    showMS1_.setVisible(false);
    addItemPopup_.add(showMS1_);

    item = new JMenuItem(LABEL_SORT_BY_ORDER);
    item.addActionListener(this);
    item.setVisible(orderType!=ORDER_TYPE_AS_IS);
    addItemPopup_.add(item);
    
    item = new JMenuItem(LABEL_SORT_BY_MASS);
    item.addActionListener(this);
    item.setVisible(orderType!=ORDER_TYPE_MZ);
    addItemPopup_.add(item);
    
    item = new JMenuItem(LABEL_SORT_BY_INTENSITY);
    item.addActionListener(this);
    item.setVisible(orderType!=ORDER_TYPE_INTENSITY);
    addItemPopup_.add(item);
    
    //Edit MSn rule
    item = new JMenuItem(LABEL_EDIT_RULES);
    item.addActionListener(this);
    addItemPopup_.add(item);

    //Recalculate an MSn identification
    if (!model.isShowMSn() && !model.isShowOmega()){
      item = new JMenuItem(LABEL_RECALCULATE_MSN);
      item.addActionListener(this);
      addItemPopup_.add(item);
      
      item = new JMenuItem(LABEL_EDIT_RT);
      item.addActionListener(this);
      addItemPopup_.add(item);
    }
    
    if (Settings.SHOW_OMEGA_TOOLS) {
      item = new JMenuItem(LABEL_EDIT_OMEGA_ASSIGNMENT);
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
              showMs2OrMs1Selection();
              addItemPopup_.show(e.getComponent(), e.getX(), e.getY());
            }
            lastSelectedIndex_ = rowNumber;
          }else{
            controlDown_ = false;
            shiftDown_ = false;
            parentListener_.listSelectionChanged(rowNumber);
            if (SwingUtilities.isRightMouseButton( e )){
              model.setSelectionInterval( rowNumber, rowNumber );
              showMs2OrMs1Selection();
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
    this.hasOh_ = model.hasOhInfo();
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
    int position;
    LipidParameterSet params;
    String displayString;
    
    switch (actionCommand) {
      
      case LABEL_ADD_ANALYTE_BEFORE:
        position = getSelectionModel().getLeadSelectionIndex();
        params = parentListener_.getAnalyteInTableAtPosition(position);
        displayString = params.getNameString();
        if (params.getModificationName()!=null&&params.getModificationName().length()>0)displayString+="_"+params.getModificationName();
        new AddAnalyteDialog(new JFrame(),"Enter new analyte", "Add a new analyte before "+displayString,params.getNameString(),params.Mz[0],params.getAnalyteFormula(),
            params.getModificationName(),params.getModificationFormula(),position,showRt_,hasOh_ ? String.valueOf(params.getOhNumber()) : null,parentListener_);
        break;
        
      case LABEL_ADD_ANALYTE_AFTER:
        position = getSelectionModel().getLeadSelectionIndex()+1;
        params = parentListener_.getAnalyteInTableAtPosition(position-1);
        displayString = params.getNameString();
        if (params.getModificationName()!=null&&params.getModificationName().length()>0)displayString+="_"+params.getModificationName();
        new AddAnalyteDialog(new JFrame(),"Enter new analyte", "Add a new analyte after "+displayString,params.getNameString(),params.Mz[0],params.getAnalyteFormula(),
            params.getModificationName(),params.getModificationFormula(),position,showRt_,hasOh_ ? String.valueOf(params.getOhNumber()) : null,parentListener_);
        break;
        
      case LABEL_DELETE_ANALYTE:
        int[] indices = getSelectedRows();
        int count = 0;
        displayString = "\n";
        for (int ind : indices){
          count++;
          params = parentListener_.getAnalyteInTableAtPosition(ind);
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
        break;
        
      case LABEL_DELETE_MOLECULAR_SPECIES:
      	ArrayList<Integer> msnIndices = new ArrayList<Integer>();
        displayString = "\n";
        for (int ind : getSelectedRows())
        {
        	params = parentListener_.getAnalyteInTableAtPosition(ind);
        	if (params instanceof LipidomicsMSnSet && (((LipidomicsMSnSet)params).getStatus()>LipidomicsMSnSet.HEAD_GROUP_DETECTED))
        	{
        		msnIndices.add(ind);
          	displayString += parentListener_.getNameInTableAtPosition(ind);
          	displayString += ", ";
            if (msnIndices.size()%8 == 0) displayString += "\n";
        	}
        }
        if (msnIndices.isEmpty())
        {
        	new WarningMessage(new JFrame(), "Error trying to delete molecular species.", 
        			String.format("The selected analytes do not contain any identifications at the molecular species level. Use the option '%s' instead!", LABEL_DELETE_ANALYTE));
        }
        else
        {
        	if (displayString.endsWith("\n")) displayString = displayString.substring(0,displayString.length()-1);
          if (displayString.length()>1) displayString = displayString.substring(0,displayString.length()-2);
          if (JOptionPane.showConfirmDialog(this, "Do you really want to delete "+displayString+"?") == JOptionPane.YES_OPTION){
            parentListener_.removeMolecularSpecies(msnIndices);
          }
        }
      	break;
        
      case LABEL_SHOW_MSMS:
        position = getSelectionModel().getLeadSelectionIndex();
        parentListener_.showMs2(position);
        break;
        
      case LABEL_SHOW_MS1:
        position = getSelectionModel().getLeadSelectionIndex();
        parentListener_.initANewViewer(position);
        break;
        
      case LABEL_SORT_BY_ORDER:
        parentListener_.changeListSorting(ORDER_TYPE_AS_IS);  
        break;
        
      case LABEL_SORT_BY_MASS:
        parentListener_.changeListSorting(ORDER_TYPE_MZ);
        break;
        
      case LABEL_SORT_BY_INTENSITY:
        parentListener_.changeListSorting(ORDER_TYPE_INTENSITY);
        break;
        
      case LABEL_EDIT_RULES:
        position = getSelectionModel().getLeadSelectionIndex();
        parentListener_.newRule(position);   
        break;
        
      case LABEL_RECALCULATE_MSN:
        position = getSelectionModel().getLeadSelectionIndex();
        parentListener_.recalculateMSn(position); 
        break;
        
      case LABEL_EDIT_RT:
        position = getSelectionModel().getLeadSelectionIndex();
        parentListener_.editRt(position);
        break;
        
      case LABEL_EDIT_OMEGA_ASSIGNMENT:
        position = getSelectionModel().getLeadSelectionIndex();
        parentListener_.editOmegaAssignment(position);
        break;
        
    }
  }
  
  public String getDisplayedNameAt(int rowIndex)
  {
  	return ((LipidomicsTableModel)super.getModel()).getDisplayedNameAt(rowIndex);
  }
  
  
//  public String getSumLipidNameAt(int rowIndex){
//    return ((LipidomicsTableModel)super.getModel()).getSumLipidNameAt(rowIndex);
//  }
  
  /**
   * depending on the currently displayed view, activates or deactivates the corresponding buttons
   * for changing the view
   */
  private void showMs2OrMs1Selection(){
    if (parentListener_.isMS2Showing()){
      this.showMS2_.setVisible(false);
      this.showMS1_.setVisible(true);      
    }else{
      this.showMS2_.setVisible(true);
      this.showMS1_.setVisible(false);
    }
      
  }
  
}
