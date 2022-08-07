package node;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.Stack;

import edge.Edge;

public class Node {
	Node root; //reference of root node
	String name; // name
	String realName; // real name for gene nodes
	String internalNodeName;//internal node's name for simulated data
	int numberOfChild; // the number of children
	LinkedList<Node>childrenList; //
	Node[] children;
	Node left;
	Node right;
	int postOrderNumber;
	int inOrderNumber;
	int preOrderNumber;
	int height = -1;
	int depth;
	Node parent;
	Node leftMostDescendent;
	Node rightMostDescendent;
	double rightBranchLength;
	double leftBranchLength;
	double parentBranchLength;
	Edge leftBranch;
	Edge rightBranch;
	String ID;
	int mappingCount=0;
	int duplicationCount=0;
	int transferCount=0;
	int transferRecipientCount=0;
	ArrayList<Node> candidateNodes = new ArrayList<Node>();
	public Node(String name){
		this.name = name;
		children = new Node[2];
	}
	public Node(Node node) {
		this.name = node.name;
		this.left = node.left;
		this.right = node.right;
		this.parent = node.parent;
	}
	
	public int getDuplicationCount() {
		return duplicationCount;
	}
	public void setDuplicationCount(int duplicationCount) {
		this.duplicationCount = duplicationCount;
	}
	public int getTransferCount() {
		return transferCount;
	}
	public void setTransferCount(int transferCount) {
		this.transferCount = transferCount;
	}
	public int getTransferRecipientCount() {
		return transferRecipientCount;
	}
	public void setTransferRecipientCount(int transferRecipientCount) {
		this.transferRecipientCount = transferRecipientCount;
	}
	public String getID() {
		return ID;
	}
	public void setID(String iD) {
		ID = iD;
	}
	public Edge getLeftBranch() {
		return leftBranch;
	}
	public void setLeftBranch(Edge leftBranch) {
		this.leftBranch = leftBranch;
	}
	public Edge getRightBranch() {
		return rightBranch;
	}
	public void setRightBranch(Edge rightBranch) {
		this.rightBranch = rightBranch;
	}
	public Node getRoot() {
		return root;
	}
	public void setRoot(Node root) {
		this.root = root;
	}
	public Node[] getChildren() {
		return children;
	}
	public void setChildren(Node[] children) {
		this.children = children;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getRealName() {
		return realName;
	}
	public void setRealName(String realName) {
		this.realName = realName;
	}
	public int getNumberOfChild() {
		return numberOfChild;
	}
	public void setNumberOfChild(int numberOfChild) {
		this.numberOfChild = numberOfChild;
	}
	public LinkedList<Node> getChildrenList() {
		return childrenList;
	}
	public void setChildrenList(LinkedList<Node> childrenList) {
		this.childrenList = childrenList;
	}
	public Node getLeft() {
		return left;
	}
	public void setLeft(Node left) {
		left.setParent(this);
		this.left = left;
		children[0]= left;
	}
	public Node getRight() {
		return right;
	}
	public void setRight(Node right) {
		right.setParent(this);
		this.right = right;
		children[1]= right;
	}
	public int getPostOrderNumber() {
		return postOrderNumber;
	}
	public void setPostOrderNumber(int postOrderNumber) {
		this.postOrderNumber = postOrderNumber;
	}
	public int getInOrderNumber() {
		return inOrderNumber;
	}
	public void setInOrderNumber(int inOrderNumber) {
		this.inOrderNumber = inOrderNumber;
	}
	public int getPreOrderNumber() {
		return preOrderNumber;
	}
	public void setPreOrderNumber(int preOrderNumber) {
		this.preOrderNumber = preOrderNumber;
	}
	public int getDepth() {
		return depth;
	}
	public void setDepth(int height) {
		this.depth = height;
	}
	public Node getParent() {
		return parent;
	}
	public void setParent(Node parent) {
		this.parent = parent;
	}
	
	public String getInternalNodeName() {
		return internalNodeName;
	}
	public void setInternalNodeName(String internalNodeName) {
		this.internalNodeName = internalNodeName;
	}
	public void setHeight(int height) {
		this.height = height;
	}
	public int getHeight() {
		if(height>-1)return height;
		if(left==null && right ==null) {
			height = 0;
		}
		else if (left ==null)height = right.getHeight()+1;
		else if (right == null)height= left.getHeight() + 1;
		else height= Math.max(left.getHeight(), right.getHeight())+1;
		return height;
	}
	public Node getLeftMostDescendent() {
		if(leftMostDescendent==null) {
			Node leftNode = left;
			
			while(leftNode.getLeft()!=null) {
				leftNode= leftNode.getLeft();
			}
			leftMostDescendent = leftNode;
		}
		return leftMostDescendent;
	}
	public void setLeftMostDescendent(Node leftMostDescendent) {
		this.leftMostDescendent = leftMostDescendent;
	}
	public Node getRightMostDescendent() {
		if(rightMostDescendent==null) {
			Node rightNode = right;
			
			while(rightNode.getRight()!=null) {
				rightNode= rightNode.getRight();
			}
			rightMostDescendent = rightNode;
		}
		return rightMostDescendent;
	}
	public void setRightMostDescendent(Node rightMostDescendent) {
		this.rightMostDescendent = rightMostDescendent;
	}
	public double getRightBranchLength() {
		return rightBranchLength;
	}
	public void setRightBranchLength(double rightBranchLength) {
		this.rightBranchLength = rightBranchLength;
	}
	public double getLeftBranchLength() {
		return leftBranchLength;
	}
	public void setLeftBranchLength(double leftBranchLength) {
		this.leftBranchLength = leftBranchLength;
	}
	public double getParentBranchLength() {
		return parentBranchLength;
	}
	public void setParentBranchLength(double parentBranchLength) {
		this.parentBranchLength = parentBranchLength;
	}

	public int getMappingCount() {
		return mappingCount;
	}
	public void setMappingCount(int mappingCount) {
		this.mappingCount = mappingCount;
	}
	public void decrementMappingCount() {
		mappingCount--;
	}
	public Node getCalculatedRightMostDescendent() {
		Node calculated = this;
		//go to right node as long as possible 
		while(calculated.getRight()!=null) {
		
			calculated = calculated.getRight();
		}
		rightMostDescendent = calculated;
		return calculated;
	}
	public Node getCalculatedLeftMostDescendent() {
		Node calculated = this;
		//go to left node as long as possible 
		while(calculated.getLeft()!=null) {
			calculated = calculated.getLeft();
		}
		leftMostDescendent = calculated;
		return calculated;
	}
	public Node getSibling() {
		Node parent = getParent();
		//root node
		if(parent==null)return null;
		//if it is left child of the parent then it should return the right child of the parent
		if(parent.getLeft()!= null && parent.getLeft().equals(this))return parent.getRight();
		//else return the left child of the parent	
		else return parent.getLeft();
	}
	//check if this node is the root of the tree;
	public boolean isRoot() {
		if(this.parent==null)return true;
		else return false;
	}
	//check if this node is the leaf of the tree;
	public boolean isLeaf() {
		//leaf node has no valid child
		if(this.left==null && this.right==null)return true;
		else return false;
	}
	//swap two children of a node
	public void swapChildren() {
		Node temp = new Node(left);
		this.left = this.right;
		this.right = temp;
	}
	//connect this node to a parent node
	public void connectAsChild(Node parent) {
		this.parent = parent;
		if(parent.left==null)parent.left = this;
		else parent.right = this;
	}
	//get list of all descendants in PreOrder of this node
	public ArrayList<Node> getDescendents(){
		ArrayList<Node> list = new ArrayList<Node>();
		Stack<Node>stack = new Stack<Node>();
		
		if(this.left!=null)stack.push(this.left);
		if(this.right!=null)stack.push(this.right);
		
		while(!stack.isEmpty()) {
			Node node = stack.pop();
			list.add(node);
			
			if(node.getRight()!=null)
			{
				stack.push(node.getRight());
			}
			if(node.getLeft()!=null)
			{
				stack.push(node.getLeft());
			}
		}
		return list;
	}
	public boolean isDescendant(Node other) {
		if(this.preOrderNumber>other.preOrderNumber && this.postOrderNumber<other.postOrderNumber)
		return true;
		else return false;
	}
	public boolean isDescendantOrEqual(Node other) {
	
		if(this.equals(other))return true;
		else if(this.preOrderNumber>other.preOrderNumber &&  this.postOrderNumber<other.postOrderNumber )
		return true;
		else return false;
	}
	// get Least common ancestor of two nodes
	public Node getLCA(Node other) {
		//if this node is descendent of the other node then other node should be LCA
		if(this.isDescendant(other)) {
			return other;
		}
		//if other node is descendent of the this node then this node should be LCA
		else if(other.isDescendant(this)) {
			return this;
		}
		//if other node
		if(this.depth<other.depth) {
			Node parent = this.parent;
			while(true) {
				if(other.isDescendant(parent))return parent;
				else parent = parent.parent;
			}
		}
		else {
			Node parent = other.parent;
			while(true) {
				if(this.isDescendant(parent))return parent;
				else parent = parent.parent;
			}
		}

	}
	//get distance between two nodes
	public int getDistance(Node other) {
		if(other.equals(this))return 0;
		//if one is descendant of other then their difference of height gives distance
		if(this.isDescendant(other)) {
			return this.depth - other.depth;
		}
		else if(other.isDescendant(this)) {
			return other.depth - this.depth;
		}
		//otherwise get LCA= Least Common Ancestor first and then calculate distance via LCA
		else {
			Node LCA = this.getLCA(other);
			return this.getDistance( LCA)+other.getDistance(LCA);
		}
	}
	
	//get all ancestors of this node by going up the tree by parent upto root
	public ArrayList<Node> getAllAncestors(){
		ArrayList<Node>ancestors = new ArrayList<>();
		if(this.parent==null)return ancestors;
		Node parent = this.parent;
		while(parent!=null) {
			ancestors.add(parent);
			parent = parent.getParent();
		}
		return ancestors;
	}
	//print node
	@Override
	public String toString() {
		return name;
	}
	public ArrayList<Node> getLosses(Node parentImage) {
		
		ArrayList<Node> losses = new ArrayList<>();
		if(this.equals(parentImage))return losses;
		Node present = this.parent;
		while(!present.equals(parentImage)) {
			losses.add(present);
			present = present.getParent();
		}
		return losses;
	}
	public ArrayList<Node> getInTrasferLosses(Node parentImage) {
		
		ArrayList<Node> losses = new ArrayList<>();
		if(this.equals(parentImage))return losses;
		Node present = this.parent;
		while(!present.equals(parentImage)) {
			losses.add(present);
			present = present.getParent();
		}
		if(parentImage.getLeft().equals(this)) {
			losses.add(parentImage.getRight());
		}
		else {
			losses.add(parentImage.getLeft());
			
		}
		return losses;
	}
	public ArrayList<Edge> getLostBranches(SpeciesNodeRooted parentImage) {
		ArrayList<Edge> lostBranches = new ArrayList<>();
		if(this.equals(parentImage))return lostBranches;
		Node present = this.parent;
		while(!present.equals(parentImage)) {
			if(this.equals(present.left)) {
				lostBranches.add(present.getRightBranch());
			}
			else {
				lostBranches.add(present.getLeftBranch());
				
			}
			present = present.getParent();
		}
		return lostBranches;
	}
	public ArrayList<Edge> getInTransferLostBranches(SpeciesNodeRooted parentImage) {
		ArrayList<Edge> lostBranches = new ArrayList<>();
		if(this.equals(parentImage))return lostBranches;
		Node present = this.parent;
		while(!present.equals(parentImage)) {
			if(this.equals(present.left)) {
				lostBranches.add(present.getRightBranch());
			}
			else {
				lostBranches.add(present.getLeftBranch());
				
			}
			present = present.getParent();
		}
		if(parentImage.getLeft().equals(this)) {
			lostBranches.add(parentImage.getRightBranch());
		}
		else {
			lostBranches.add(parentImage.getLeftBranch());
			
		}
		return lostBranches;
	}
	public void setBranches() {
		if(isLeaf())return;
		leftBranch = new Edge(this, left);
		rightBranch = new Edge(this, right);
	}
	public String getLCAFormat() {
		String leftNodeName = getLeftMostDescendent().name;
		String rightNodeName = getRightMostDescendent().name;
		return "LCA["+leftNodeName+","+rightNodeName+"]";
	}
	public boolean isOnLeftBranch() {
		if(parent==null)return false;
		if(parent.left!=null && parent.getLeft().equals(this)) 
			return true;
		
		else return false;
	}
	public boolean isOnRightBranch() {
		if(parent==null)return false;
		if(parent.right!=null && parent.getRight().equals(this)) 
			return true;
		
		else return false;
	}
	public Edge getParentBranch() {
		if(parent==null)return null;
		else if(isOnLeftBranch())return parent.leftBranch;
		else return parent.rightBranch;	
	}
	public ArrayList<Node> getAllLeavesUnderSubtree() {
		ArrayList<Node> leavesUnderSubtree = new ArrayList<Node>();
		ArrayList<Node> allDecendentNodes = getDescendents();
		allDecendentNodes.add(this);
		for(Node node: allDecendentNodes) {
			if(node.isLeaf()) {
				leavesUnderSubtree.add(node);
			}
		}
		return leavesUnderSubtree;
	}
	public boolean anyAncestorBranchHavingLosses() {
		Edge parentBranch = getParentBranch();
		Node parentPointer = parent;
		int i =1;
	
		
		while(parentBranch!=null && i<3 ) {
			if(parentBranch.getNumberOfLoss()>0) {
				return true;
				
			}
			parentBranch = parentBranch.getParent().getParentBranch();
			parentPointer = parentPointer.parent;
			i++;
		}
		return false;
	}
	public ArrayList<Node> getSiblingLeaves(){
		ArrayList<Node> siblingLeaves = new ArrayList<Node>();
		Node sibling = getSibling();
		if(this.isRoot() || parent.isRoot() || sibling==null)return siblingLeaves;
		siblingLeaves = sibling.getAllLeavesUnderSubtree();
		return siblingLeaves;
	}
	public ArrayList<Node> getUncleLeaves(){
		ArrayList<Node> uncleLeaves = new ArrayList<Node>();
		if(this.isRoot() || this.parent.isRoot())return uncleLeaves;
		Node sibling = getParent().getSibling();
		uncleLeaves = sibling.getAllLeavesUnderSubtree();
		return uncleLeaves;
	}
	public ArrayList<Node> getAllKindOfSiblingLeaves() {
		ArrayList<Node> siblingLeaves  =	getSiblingLeaves();
	
		ArrayList<Node> uncleLeaves = getUncleLeaves();
		 siblingLeaves.addAll(uncleLeaves);
	
		return siblingLeaves;
	}
	public Node getNthParent(int heightOfNode) {
		Node node = this;
		int i = 0;
		while(i<heightOfNode) {
			if(node.parent==root)break;
			node= node.parent;
			i++;
		}
		return node;
	}
	public void deleteAsLeafFromTree() {
		if(isOnLeftBranch()) {
			parent.left = null;
		}
		else parent.right = null;
	}
	public void deleteAsSingleDegreeInternalNode() {
		if(isRoot())return;
		Node child = null;
		if(this.left==null)child = right;
		else child = left;
		if(isOnLeftBranch()) {
			parent.left = child;
			child.parent = parent;
		}
		else {
			parent.right = child;
			child.parent = parent;
		}
	}
	public void incrementMappingCount() {
		mappingCount++;
	}
	public boolean hasPositiveParentMapping() {
		Node parentCursor = parent;
		
		while(parentCursor!=null ) {
			if(parentCursor.mappingCount>0) {
				
				return true;
				
			}
			parentCursor = parentCursor.parent;
		
		}
		return false;
		
	}
	public void incrementDuplicationCount() {
		duplicationCount++;
	}
	public void incrementTransferCount() {
		transferCount++;
	}
	public void incrementTransferRecipientCount() {
		transferRecipientCount++;
		
	}
	public int getInferredFrequencyByCurrentHeuristic(Node rootImage) {
		int sumOfDuplicationCounts = 0;
		int sumOfTransferRecipientCounts = 0;
		int rootMappingCount = 0;
		
		if(this.isDescendantOrEqual(rootImage))rootMappingCount = 1;
		Node presentNodeCursor = this;
		
		while(presentNodeCursor!=null ) {
			sumOfDuplicationCounts+= presentNodeCursor.duplicationCount;
			sumOfTransferRecipientCounts+= presentNodeCursor.transferRecipientCount;
			presentNodeCursor = presentNodeCursor.parent;
		
		}
	//	System.out.println("Leaf "+ this.getName());
	//	System.out.println("Dup: "+sumOfDuplicationCounts);
	//	System.out.println("Rec: "+sumOfTransferRecipientCounts);
	//	System.out.println("Root: "+ rootMappingCount);
		
		int total = sumOfDuplicationCounts + sumOfTransferRecipientCounts+ rootMappingCount;
//		System.out.println("Inferred count: "+ total);
	//	System.out.println();
		return total;
	}
	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + postOrderNumber;
		return result;
	}
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Node other = (Node) obj;
		if (name == null) {
			if (other.name != null)
				return false;
		} else if (!name.equals(other.name))
			return false;
		if (postOrderNumber != other.postOrderNumber)
			return false;
		return true;
	}
	public int getTotalDuplicationCount() {
		int sumOfDuplicationCounts = 0;
		
		Node presentNodeCursor = this;
		
		while(presentNodeCursor!=null ) {
			sumOfDuplicationCounts+= presentNodeCursor.duplicationCount;
			presentNodeCursor = presentNodeCursor.parent;
		
		}
		
		
	return sumOfDuplicationCounts;
	
	}
	public int getTotalLossCount() {
		Edge parentBranch = getParentBranch();
		Node parentPointer = parent;
		int loss =0;
		while(parentBranch!=null ) {
			loss+= parentBranch.getNumberOfLoss();
			parentBranch = parentBranch.getParent().getParentBranch();
			parentPointer = parentPointer.parent;
		
		}
	
		return loss;
	}

	private boolean isDescendentOrAncestor(Node anotherNode) {
		if(this.isDescendantOrEqual(anotherNode) || anotherNode.isDescendantOrEqual(this))
		return true;
		else return false;
	}
	public void addCandidateNode(GeneNodeRooted geneNode) {
		Node gene = geneNode ;
		candidateNodes.add(gene);
	}
	
	public boolean anyAncestorBranchHavingLosses(int i) {
		Edge parentBranch = getParentBranch();
		Node parentPointer = parent;
		int level = 1;
	
		
		while(parentBranch!=null && level<=i ) {
			if(parentBranch.getNumberOfLoss()>0) {
				return true;
				
			}
			parentBranch = parentBranch.getParent().getParentBranch();
			parentPointer = parentPointer.parent;
			i++;
		}
		return false;
	}
	public boolean nonDescendentOfRootImage(ArrayList<Node> possibleRootImages) {
		for(Node rootImage: possibleRootImages)
		{
			if(rootImage.isDescendentOrAncestor(this))
				return false;
		}
		return true;
	}
	public boolean additiveByMappingCount(Node node) {
	
		for(Node candidateNode: candidateNodes) {
			if(!candidateNode.isDescendantOrEqual(node) && !node.isDescendantOrEqual(candidateNode))
				return true;
		}
		return false;
	}
	
}
