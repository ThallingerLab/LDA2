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

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.Vector;

import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.ListSelectionModel;

import at.tugraz.genome.lda.TooltipTexts;
import at.tugraz.genome.lda.utils.StaticUtils;

/**
 * 
 * @author Juergen Hartler
 *
 */
public class GroupPanel extends JPanel
{

  private static final long serialVersionUID = -7820284034078799038L;
  
  private String groupName_;
  private Vector<File> selectedFiles_ = new Vector<File>();
  private JTable selectedTable_;
  private int maxWidth_;
  ActionListener parentListener_;
  private ImageIcon editIcon_ = new ImageIcon(getClass().getResource(
  "/images/Edit.gif"));
  private ImageIcon removeFilesIcon_ = new ImageIcon(getClass().getResource(
  "/images/removeFiles.gif"));
  private ImageIcon deleteIcon_ = new ImageIcon(getClass().getResource(
  "/images/Delete.gif"));
  
  public GroupPanel(String groupName,Vector<File> selectedFiles, int maxWidth, ActionListener parentListener){
    super();
    this.groupName_ = groupName;
    setLayout(new BorderLayout());   
    selectedFiles_ = selectedFiles;
    this.maxWidth_ = maxWidth;
    this.parentListener_ = parentListener;
    this.createComponents();
  }
  
  private void createComponents(){
    this.createTopMenu();
    this.generateTable();  
  }
  
  private void createTopMenu(){
    JPanel topMenuPanel = new JPanel();
    topMenuPanel.setLayout(new GridBagLayout());
    
    JLabel name = new JLabel(groupName_);
    name.setToolTipText(TooltipTexts.STATISTICS_GROUP_SELECTION_NAME);
    topMenuPanel.add(name,new GridBagConstraints(0, 0, 3, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
    JButton rename = new JButton("Rename",editIcon_);
    rename.setIconTextGap(2);
    rename.setMargin(new Insets(2,2,2,2));
    rename.addActionListener(parentListener_);
    rename.setActionCommand("renameGroup;"+groupName_);
    rename.setToolTipText(TooltipTexts.STATISTICS_GROUP_SELECTION_RENAME);
    topMenuPanel.add(rename ,new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    JButton remove = new JButton("Remove",removeFilesIcon_);
    remove.setIconTextGap(2);
    remove.setMargin(new Insets(2,2,2,2));
    remove.addActionListener(parentListener_);
    remove.setActionCommand("removeItems;"+groupName_);
    remove.setToolTipText(TooltipTexts.STATISTICS_REMOVE_SELECTION);
    topMenuPanel.add(remove ,new GridBagConstraints(1, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    JButton delete = new JButton("Delete Group",deleteIcon_);
    delete.setIconTextGap(2);
    delete.setMargin(new Insets(2,2,2,2));
    delete.addActionListener(parentListener_);
    delete .setActionCommand("deleteGroup;"+groupName_);
    delete.setToolTipText(TooltipTexts.STATISTICS_GROUP_SELECTION_DELETE);
    topMenuPanel.add(delete ,new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

    
    this.add(topMenuPanel,BorderLayout.NORTH);
  }

  private void generateTable(){
    String[][] tableData = new String[selectedFiles_.size()][1];
//    int longestFirstColumnElement = 0;
    for (int i=0; i!=selectedFiles_.size();i++){
      String fileName = StaticUtils.extractFileName(selectedFiles_.get(i).getAbsolutePath());
      tableData[i][0] = fileName;
//      if (fileName.length()>longestFirstColumnElement)longestFirstColumnElement = fileName.length();
    }
    String[] columnNames = { "file name"};
    selectedTable_ = new JTable(tableData, columnNames);
    ListSelectionModel selectionModel = selectedTable_.getSelectionModel();
    selectedTable_.setSelectionModel(selectionModel);
    selectionModel.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);  
    JScrollPane scrollPane = new JScrollPane(selectedTable_);
    scrollPane.setPreferredSize(new Dimension(maxWidth_, 95));
    if (tableData!=null && tableData.length>0){
//      int longestSecondColumnElement = 0;
//      for (int i=0;i!=tableData.length;i++){
//        int lengthOne = tableData[i][0].length();
        
//      }
//      selectedTable_.getColumnModel().getColumn(0).setPreferredWidth((int)(this.maxWidth_*0.95));
//      selectedTable_.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);
    }
    this.add(scrollPane,BorderLayout.CENTER);
  }

  public JTable getTable(){
    return this.selectedTable_;
  }
  
  public void updateFiles(Vector<File> files){
    selectedFiles_ = files;
    this.removeAll();
    this.createComponents();
    this.invalidate();
    this.updateUI();
  }

  
}
