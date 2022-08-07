package node;

import mapping.RootedBinaryTreeMapping;

public class GeneNodeRooted extends Node{
	private SpeciesNodeRooted mappedSpeciesNode;
	private int option;
	private String event;
	private RootedBinaryTreeMapping mapping;
	public GeneNodeRooted(Node node) {
		super(node);
	}
	public GeneNodeRooted(String name){
		super(name);
	}
	public SpeciesNodeRooted getMappedSpeciesNode() {
		return mappedSpeciesNode;
	}
	public void setMappedSpeciesNode(SpeciesNodeRooted mappedSpeciesNode) {
		this.mappedSpeciesNode = mappedSpeciesNode;
	}
	public int getOption() {
		return option;
	}
	public void setOption(int option) {
		this.option = option;
	}
	public String getEvent() {
		return event;
	}
	public void setEvent(String event) {
		this.event = event;
	}
	
	public RootedBinaryTreeMapping getMapping() {
		return mapping;
	}
	public void setMapping(RootedBinaryTreeMapping mapping) {
		this.mapping = mapping;
	}
	public void setEventByOption() {
		switch(option) {
		case 0: 
			event = "Leaf Node";
			break;
		case 1: 
			event = "Duplication";
			break;
		case 2: 
		
		case 3: 
			event = "Transfer";
			break;
		case 4: 
			
		case 5: 
			event = "Speciation";
			break;
		case 6: 
			
		case 7: 
		case 8:
		case 9:
		case 10:
		case 11:
		case 12:
		case 13:
		case 14:
			event = "Duplication";
			break;
		case 15: 
			
		case 16: 
			event = "Transfer";
			break;
			
		}
		
	}
	public String getFullName() {
		if(realName!=null)
		return name+"_"+realName;
		else return name;
	}
	public void setInternalNodeName(String gene) {
		super.setInternalNodeName(gene);
	}
	
}
