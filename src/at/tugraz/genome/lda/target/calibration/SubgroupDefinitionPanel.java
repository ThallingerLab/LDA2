package at.tugraz.genome.lda.target.calibration;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Collections;

import javax.swing.AbstractCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellEditor;
import javax.swing.table.TableCellRenderer;

import at.tugraz.genome.lda.WarningMessage;
import at.tugraz.genome.lda.target.JOptionPanel;
import at.tugraz.genome.lda.verifier.StringVerifier;


public class SubgroupDefinitionPanel extends JFrame
{
	private static final long serialVersionUID = 1L;
	private static final String TITLE = "Define subgroups to be calibrated together";
	
	private ArrayList<SubGroup> definedSubgroups_;
	private ArrayList<String> usedLipidClasses_;
	private ArrayList<String> ungroupedLipidClasses_;
	private DefaultTableModel defineNewTableModel_;
	private JTextField groupNameField_;
	private JPanel definedGroupsPanel_;
	private JPanel defineNewGroupPanel_;
	private JPanel contentsPanel_;
	
	
	public SubgroupDefinitionPanel(ArrayList<SubGroup> definedSubgroups, String[] displayedLipidClasses, JPanel parent)
	{
		super(TITLE);
		this.definedSubgroups_ = definedSubgroups;
		this.usedLipidClasses_ = new ArrayList<String>();
		for (SubGroup subGroup : definedSubgroups_)
		{
			usedLipidClasses_.addAll(subGroup.getLipidClasses());
		}
		this.ungroupedLipidClasses_  = new ArrayList<String>();
		for (int i=0;i<displayedLipidClasses.length;i++)
		{
			if (!usedLipidClasses_.contains(displayedLipidClasses[i]) && !displayedLipidClasses[i].equals(CalibrationGraphPanel.PLOT_ALL))
			{
				ungroupedLipidClasses_.add(displayedLipidClasses[i]);
			}
		}
		this.contentsPanel_ = new JPanel();
		contentsPanel_.setLayout(new GridBagLayout());
		this.add(contentsPanel_);
		
		initDefinedGroupsPanel();
		initDefineNewGroupPanel();
		
		//general settings
		setLocationRelativeTo(parent);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    setSize(950, 850);
    setVisible(true);
		
	}
	
	private void initDefinedGroupsPanel()
	{
		definedGroupsPanel_ = new JPanel();
		definedGroupsPanel_.setLayout(new GridBagLayout());
		
		if (definedSubgroups_.isEmpty())
		{
			JLabel label = new JLabel("No subgroups have been defined.");
			label.setPreferredSize(new Dimension(663,50));
			label.setMinimumSize(label.getPreferredSize());	
			definedGroupsPanel_.add(label, new GridBagConstraints(0, 0, 0, 1, 0.0, 0.0
	        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(25, 25, 25, 25), 0, 0));
		}
		else
		{
			String[] headers = { "Defined subgroups ", "Actions" };
			
			DefaultTableModel model = new DefaultTableModel()
			{
	      private static final long serialVersionUID = 1L;
	      
	      @Override
	      public boolean isCellEditable(int row, int column) 
	      {
	      	switch (column) 
	      	{
	          case 1:
	            return true;
	          default:
	            return false;
	      	}
	      }
			};
			for (String header : headers)
			{
				model.addColumn(header);
			}
			for (SubGroup definedSubgroup : definedSubgroups_) 
	    {
				Object[] row = { definedSubgroup.getId(), new RemoveAction(false, definedSubgroup, this) };
				model.addRow(row);
	    }
			JTable table = new JTable(model);
			table.getColumnModel().getColumn(1).setCellEditor(new RemoveActionEditor());
			table.getColumnModel().getColumn(1).setCellRenderer(new RemoveActionRenderer());
			
			JScrollPane pane = new JScrollPane(table);
			pane.setPreferredSize(new Dimension(700,125));
			pane.setMinimumSize(pane.getPreferredSize());
			definedGroupsPanel_.add(pane, new GridBagConstraints(0, 0, 0, 1, 0.0, 0.0
	        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		}
		
		definedGroupsPanel_.setBorder(JOptionPanel.getTitledPanelBorder("Current subgroups"));
		
		contentsPanel_.add(definedGroupsPanel_, new GridBagConstraints(0, 0, 0, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(10, 25, 0, 25), 0, 0));
	}
	
	private void removeSubgroup(SubGroup subGroup)
	{
		definedSubgroups_.remove(subGroup);
		usedLipidClasses_.removeAll(subGroup.getLipidClasses());
		ungroupedLipidClasses_.addAll(subGroup.getLipidClasses());
		contentsPanel_.remove(definedGroupsPanel_);
		contentsPanel_.remove(defineNewGroupPanel_);
		initDefinedGroupsPanel();
		initDefineNewGroupPanel();
		contentsPanel_.invalidate();
		contentsPanel_.updateUI();
	}
	
	private void initDefineNewGroupPanel()
	{
		defineNewGroupPanel_ = new JPanel();
		defineNewGroupPanel_.setLayout(new GridBagLayout());
		
		initDefineNewTableModel();
		JScrollPane pane = new JScrollPane(new JTable(defineNewTableModel_));
		pane.setPreferredSize(new Dimension(450,300));
		pane.setMinimumSize(pane.getPreferredSize());
		defineNewGroupPanel_.add(pane, new GridBagConstraints(0, 0, 1, 3, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		
		JSeparator verticalSeparator = new JSeparator(SwingConstants.VERTICAL);
		verticalSeparator.setPreferredSize(new Dimension(10,(int) pane.getPreferredSize().getHeight()));
		verticalSeparator.setMinimumSize(verticalSeparator.getPreferredSize());
		verticalSeparator.setForeground(new Color(105,123,140));
		defineNewGroupPanel_.add(verticalSeparator, new GridBagConstraints(1, 0, 1, 3, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		
		JLabel label = new JLabel("Group name: ");
		defineNewGroupPanel_.add(label, new GridBagConstraints(2, 0, 1, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		
		groupNameField_ = new JTextField(20);
		groupNameField_.setInputVerifier(new StringVerifier(1,20,"Group name: "));
		groupNameField_.setMinimumSize(groupNameField_.getPreferredSize());
		defineNewGroupPanel_.add(groupNameField_, new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
        ,GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		
		JButton button = new JButton("Save Subgroup");
		button.addActionListener(new ActionListener() {
		  public void actionPerformed(ActionEvent e) 
		  {
		  	saveSubgroup_actionPerformed(e);
		  }
	  });
		defineNewGroupPanel_.add(button, new GridBagConstraints(2, 3, 1, 1, 0.0, 0.0
        ,GridBagConstraints.SOUTHEAST, GridBagConstraints.NONE, new Insets(5, 5, 5, 5), 0, 0));
		
		defineNewGroupPanel_.setBorder(JOptionPanel.getTitledPanelBorder("Define new subgroup"));
		contentsPanel_.add(defineNewGroupPanel_, new GridBagConstraints(0, 1, 0, 1, 0.0, 0.0
        ,GridBagConstraints.CENTER, GridBagConstraints.NONE, new Insets(25, 25, 25, 25), 0, 0));
	}
	
	private void initDefineNewTableModel()
	{
		String[] headers = { "Ungrouped lipid classes", "Selected" };
		defineNewTableModel_ = new DefaultTableModel()
		{
      private static final long serialVersionUID = 1L;
      
      @Override
      public boolean isCellEditable(int row, int column) 
      {
      	switch (column) 
      	{
          case 1:
            return true;
          default:
            return false;
      	}
      }
      
      @Override
      public Class<?> getColumnClass(int column) 
      {
        switch (column) 
        {
          case 1:
            return Boolean.class;
          default:
            return String.class;
        }
      }
    };
		for (String header : headers)
		{
			defineNewTableModel_.addColumn(header);
		}
		for (String ungrouped : ungroupedLipidClasses_) 
    {
			Object[] row = { ungrouped, false };
			defineNewTableModel_.addRow(row);
    }
	}
	
	private void saveSubgroup_actionPerformed(ActionEvent e)
  {
		ArrayList<String> lipidClasses = new ArrayList<String>();
		String groupName = groupNameField_.getText();
		
		for (int i=0; i<defineNewTableModel_.getRowCount(); i++)
		{
			if ((boolean) defineNewTableModel_.getValueAt(i, 1) == true)
			{
				lipidClasses.add((String) defineNewTableModel_.getValueAt(i, 0));
			}
		}
		if (groupName.length() > 0 && lipidClasses.size()>1)
		{
			definedSubgroups_.add(new SubGroup(groupName, lipidClasses));
			usedLipidClasses_.addAll(lipidClasses);
			ungroupedLipidClasses_.removeAll(lipidClasses);
			
			contentsPanel_.remove(definedGroupsPanel_);
			contentsPanel_.remove(defineNewGroupPanel_);
			initDefinedGroupsPanel();
			initDefineNewGroupPanel();
			contentsPanel_.invalidate();
			contentsPanel_.updateUI();
		}
		else
		{
			new WarningMessage(new JFrame(), "Error", "A group name must be defined and at least 2 lipid classes must be selected!");
		}
  }
	
	public ArrayList<SubGroup> getDefinedSubgroups()
	{
		return definedSubgroups_;
	}

	public ArrayList<String> getUsedLipidClasses()
	{
		return usedLipidClasses_;
	}

	public ArrayList<String> getUngroupedLipidClasses()
	{
		Collections.sort(ungroupedLipidClasses_);
		return ungroupedLipidClasses_;
	}






	private class RemoveActionRenderer extends JCheckBox implements TableCellRenderer 
	{
		private static final long serialVersionUID = 1L;

		@Override
	  public Component getTableCellRendererComponent(
	    JTable table, Object value, boolean isSelected,
	    boolean hasFocus, int row, int col) {
			RemoveAction action = (RemoveAction) value;
			this.setSelected(action.getSelected());
	    this.setText("Remove");
	    if (action.getSelected())
	    {
	    	DefaultTableModel model = (DefaultTableModel)table.getModel();
	    	model.removeRow(row);
	    	action.getParent().removeSubgroup(action.getSubGroup());
	    }
	    return this;
	  }
	}

	private class RemoveActionEditor extends AbstractCellEditor implements TableCellEditor, ItemListener 
	{
		private static final long serialVersionUID = 1L;
		private RemoveActionRenderer renderer_ = new RemoveActionRenderer();
	  private RemoveAction action_ = null;
	
	  public RemoveActionEditor() 
	  {
	  	renderer_.addItemListener(this);
	  }
	
	  @Override
	  public Object getCellEditorValue() 
	  {
	  	action_.setSelected(renderer_.isSelected());
	    return action_;
	  }
	
	  @Override
	  public Component getTableCellEditorComponent(JTable table,
	      Object value, boolean isSelected, int row, int col) 
	  {
			action_ = (RemoveAction) value;
			renderer_.setSelected(action_.getSelected());
			renderer_.setText("Remove");
	    return renderer_;
	  }
	
	  @Override
	  public void itemStateChanged(ItemEvent e) 
	  {
	    this.fireEditingStopped();
	  }
	}

	private class RemoveAction 
	{
	  private Boolean selected_;
	  private SubGroup subGroup_;
	  private SubgroupDefinitionPanel parent_;
	
	  public RemoveAction(Boolean selected, SubGroup subGroup, SubgroupDefinitionPanel parent) 
	  {
	    this.selected_ = selected;
	    this.subGroup_ = subGroup;
	    this.parent_ = parent;
	  }
	
		public Boolean getSelected()
		{
			return selected_;
		}
	
		public void setSelected(Boolean selected)
		{
			this.selected_ = selected;
		}
	
		public SubGroup getSubGroup()
		{
			return subGroup_;
		}
	  
	  public SubgroupDefinitionPanel getParent()
	  {
	  	return this.parent_;
	  }
	}
}


