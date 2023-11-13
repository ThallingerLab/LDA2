package at.tugraz.genome.lda.target.calibration;

import java.util.ArrayList;

public class SubGroup
{
	private String groupName_;
	private ArrayList<String> lipidClasses_;
	
	public SubGroup(String groupName, ArrayList<String> lipidClasses)
	{
		this.groupName_ = groupName;
		this.lipidClasses_ = lipidClasses;
	}

	public String getGroupName()
	{
		return groupName_;
	}

	public ArrayList<String> getLipidClasses()
	{
		return lipidClasses_;
	}
	
	public String getId()
	{
		StringBuilder lipidClassesBuilder = new StringBuilder();
		for (int i=0;i<lipidClasses_.size();i++)
		{
			lipidClassesBuilder.append(lipidClasses_.get(i));
			if (i+1<lipidClasses_.size())
			{
				lipidClassesBuilder.append(", ");
			}
		}
		return String.format("%s; lipid classes %s", groupName_, lipidClassesBuilder.toString()); 
	}
}
