package mapping;

import java.math.BigInteger;
import java.util.*;

import dtl_algo.DTLReconciliation;
import edge.Edge;
import node.GeneNodeRooted;
import node.Node;
import node.SpeciesNodeRooted;

public class RootedBinaryTreeMapping {
	GeneNodeRooted geneNode;
	SpeciesNodeRooted speciesNode;
	BigInteger numberOfSolutions = new BigInteger("0");
	ArrayList<Node> inSpeciesNodes= new ArrayList<Node>();
	ArrayList<Node> outSpeciesNodes= new ArrayList<Node>();
	ArrayList<Node> inAltSpeciesNodes= new ArrayList<Node>();
	int cost = Integer.MAX_VALUE / 10;
	int inCost = Integer.MAX_VALUE / 10;
	int inAltCost = Integer.MAX_VALUE / 10;
	int outCost = Integer.MAX_VALUE / 10;
	int speciationCost = Integer.MAX_VALUE / 10;
	int duplicationCost = Integer.MAX_VALUE / 10;
	int transferCost = Integer.MAX_VALUE / 10;
	ArrayList<Integer> options = new ArrayList<Integer>();
	ArrayList<BigInteger>optionWeights = new ArrayList<>();
	ArrayList<Edge> lostEdges = new ArrayList<>();
	ArrayList<Node> losses = new ArrayList<Node>();
	Node additiveTransferRecipient = null;
	
	public RootedBinaryTreeMapping(Node geneNode, Node speciesNode) {
		this.geneNode = (GeneNodeRooted) geneNode;
		this.speciesNode = (SpeciesNodeRooted) speciesNode;
	}

	public BigInteger getNumberOfSolutions() {
		return numberOfSolutions;
	}

	public void setNumberOfSolutions(BigInteger numberOfSolutions) {
		this.numberOfSolutions = numberOfSolutions;
	}

	

	public ArrayList<Edge> getLostEdges() {
		return lostEdges;
	}

	public void setLostEdges(ArrayList<Edge> lostEdges) {
		this.lostEdges = lostEdges;
	}

	public Node getAdditiveTransferRecipient() {
		return additiveTransferRecipient;
	}

	public void setAdditiveTransferRecipient(Node additiveTransferRecipient) {
		this.additiveTransferRecipient = additiveTransferRecipient;
	}

	



	public ArrayList<Node> getOutSpeciesNodes() {
		return outSpeciesNodes;
	}

	public void setOutSpeciesNodes(ArrayList<Node> outSpeciesNodes) {
		this.outSpeciesNodes = outSpeciesNodes;
	}

	public ArrayList<Node> getInAltSpeciesNodes() {
		return inAltSpeciesNodes;
	}

	public void setInAltSpeciesNodes(ArrayList<Node> inAltSpeciesNodes) {
		this.inAltSpeciesNodes = inAltSpeciesNodes;
	}

	public ArrayList<Node> getInSpeciesNodes() {
		return inSpeciesNodes;
	}

	public void setInSpeciesNodes(ArrayList<Node> inSpeciesNodes) {
		this.inSpeciesNodes = inSpeciesNodes;
	}

	public ArrayList<Integer> getOptions() {
		return options;
	}

	public void setOptions(ArrayList<Integer> options) {
		this.options = options;
	}

	

	public GeneNodeRooted getGeneNode() {
		return geneNode;
	}

	public void setGeneNode(GeneNodeRooted geneNode) {
		this.geneNode = geneNode;
	}

	public SpeciesNodeRooted getSpeciesNode() {
		return speciesNode;
	}

	public void setSpeciesNode(SpeciesNodeRooted speciesNode) {
		this.speciesNode = speciesNode;
	}

	public ArrayList<BigInteger> getOptionWeights() {
		return optionWeights;
	}

	public void setOptionWeights(ArrayList<BigInteger> optionWeights) {
		this.optionWeights = optionWeights;
	}

	public int getCost() {
		return cost;
	}

	public void setCost(int cost) {
		this.cost = cost;
	}

	public int getInCost() {
		return inCost;
	}

	public void setInCost(int inCost) {
		this.inCost = inCost;
	}

	public int getInAltCost() {
		return inAltCost;
	}

	public void setInAltCost(int inAltCost) {
		this.inAltCost = inAltCost;
	}

	public int getOutCost() {
		return outCost;
	}

	public void setOutCost(int outCost) {
		this.outCost = outCost;
	}

	public int getSpeciationCost() {
		return speciationCost;
	}

	public void setSpeciationCost(int speciationCost) {
		this.speciationCost = speciationCost;
	}

	public int getDuplicationCost() {
		return duplicationCost;
	}

	public void setDuplicationCost(int duplicationCost) {
		this.duplicationCost = duplicationCost;
	}

	public int getTransferCost() {
		return transferCost;
	}

	public void setTransferCost(int transferCost) {
		this.transferCost = transferCost;
	}

	public ArrayList<Node> getLosses() {
		return losses;
	}

	public void setLosses(ArrayList<Node> losses) {
		this.losses = losses;
	}



	



	public void updateLosses(ArrayList<Node> nodeList1,
			ArrayList<Node> nodeList2) {
		losses.addAll(nodeList1);
		losses.addAll(nodeList2);
		
	}

	public void printLossInformation() {
		System.out.print("The number of losses: "+losses.size()+" and the lost nodes are: ");
		for(Node node: losses) {
			System.out.print(node.getName()+" ");
		}
		System.out.println();
		System.out.print("The number of lost branches: "+lostEdges.size()+" and the lost branches are: ");
		for(Edge edge: lostEdges) {
			System.out.print(edge+" ");
		}
		System.out.println();
	}

	@Override
	public String toString() {
		return   geneNode.getName() +" " + speciesNode.getName();
	}

	
	public void updateLossesWithBranch(Node childImage) {
		losses.addAll (childImage.getLosses(speciesNode));
		concatenateBranchLosses( childImage.getLostBranches(speciesNode));
		
	}
	public void addLoss(Node node) {
		losses.add(node);
		
	}

	public void updateLosses(SpeciesNodeRooted leftChildImage, SpeciesNodeRooted rightChildImage) {
		ArrayList<Node> leftChildLosses = leftChildImage.getLosses(speciesNode);
		ArrayList<Node>rightChildLosses = rightChildImage.getLosses(speciesNode);
		
		concatenateLosses( leftChildLosses);
		
		concatenateLosses( rightChildLosses);

	}

	public void updateLossesWithBranch(SpeciesNodeRooted leftChildImage, SpeciesNodeRooted rightChildImage) {
		ArrayList<Node> leftChildLosses = leftChildImage.getLosses(speciesNode);
		ArrayList<Edge> leftChildLostBranches = leftChildImage.getLostBranches(speciesNode);
		ArrayList<Node>rightChildLosses = rightChildImage.getLosses(speciesNode);
		ArrayList<Edge> rightChildLostBranches = rightChildImage.getLostBranches(speciesNode);
		
		concatenateLosses( leftChildLosses);
		
		concatenateLosses( rightChildLosses);

		concatenateBranchLosses(leftChildLostBranches);
		concatenateBranchLosses(rightChildLostBranches);
		
	}
	
	private void concatenateLosses( ArrayList<Node> newLosses) {
		for(Node node: newLosses) {
			losses.add(node);
		}
	}
	private void concatenateBranchLosses( ArrayList<Edge> newLosses) {
		for(Edge edge: newLosses) {
		addLostBranch(edge);
		
		}
	}
	public void addLostBranch(Edge branch) {
		lostEdges.add(branch);
		branch.incrementLossCount();
	}

	public void updateInTransferLossesWithBranch(SpeciesNodeRooted childImage) {
		losses.addAll(childImage.getInTrasferLosses(speciesNode));
		
		concatenateBranchLosses(childImage.getInTransferLostBranches(speciesNode));
		
	}

	public Node getInSpeciesNodeRandomly() {
		if(inSpeciesNodes.size()==0)return null;
		Random random = new Random();
		int index = random.nextInt(inSpeciesNodes.size());
		return inSpeciesNodes.get(index);
	}

	public void addInSpeciesNode(Node node) {
		inSpeciesNodes.add(node);
	}

	public void addInAltSpeciesNode(Node speciesNode) {
		inAltSpeciesNodes.add(speciesNode);
	}

	public Node getInAltSpeciesNodeRandomly() {
		if(inAltSpeciesNodes.size()==0)return null;
		Random random = new Random();
		int index = random.nextInt(inAltSpeciesNodes.size());
		return inAltSpeciesNodes.get(index);
	}

	public Node getOutSpeciesNodeRandomly() {
		if(outSpeciesNodes.size()==0)return null;
		Random random = new Random();
		int index = random.nextInt(outSpeciesNodes.size());
		return outSpeciesNodes.get(index);
	}

	public void addOutSpeciesNode(Node node) {
		if(node!=null)
		outSpeciesNodes.add(node);
		
	}

	public void addWeightByChildren(RootedBinaryTreeMapping leftToParentMapping,
			RootedBinaryTreeMapping rightToParentMapping) {
		BigInteger newWeight = leftToParentMapping.numberOfSolutions.multiply(rightToParentMapping.numberOfSolutions);
		optionWeights.add(newWeight);
		
	}

	public void addWeightByInOut(RootedBinaryTreeMapping firstMapping,
			RootedBinaryTreeMapping secondMapping, DTLReconciliation reconciliation) {
		BigInteger inCounts = new BigInteger("0");
		BigInteger outCounts = new BigInteger("0");
		for(Node inNode: firstMapping.getInSpeciesNodes()) {
			Node geneNode = firstMapping.geneNode;
			RootedBinaryTreeMapping inMapping = reconciliation.getMappingFromNodes(geneNode, inNode);
			inCounts = inCounts.add(inMapping.numberOfSolutions);
		}
		for(Node outNode: secondMapping.getOutSpeciesNodes()) {
			Node geneNode = secondMapping.geneNode;
			RootedBinaryTreeMapping outMapping = reconciliation.getMappingFromNodes(geneNode, outNode);
			outCounts = outCounts.add(outMapping.numberOfSolutions);
		}
		
		optionWeights.add(inCounts.multiply(outCounts));
		
	}

	public void addWeightByInAndIn(RootedBinaryTreeMapping firstMapping,
			RootedBinaryTreeMapping secondMapping, DTLReconciliation reconciliation) {
		BigInteger firstInCount = new BigInteger("0");
		BigInteger secondInCount = new BigInteger("0");
		for(Node inNode: firstMapping.getInSpeciesNodes()) {
			Node geneNode = firstMapping.geneNode;
			RootedBinaryTreeMapping inMapping = reconciliation.getMappingFromNodes(geneNode, inNode);
			firstInCount =firstInCount.add(inMapping.numberOfSolutions);
		}
		for(Node inNode: secondMapping.getInSpeciesNodes()) {
			Node geneNode = secondMapping.geneNode;
			RootedBinaryTreeMapping inMapping = reconciliation.getMappingFromNodes(geneNode, inNode);
			secondInCount =secondInCount.add(inMapping.numberOfSolutions);
		}
		optionWeights.add(firstInCount.multiply(secondInCount));
		
	}

	public void addWeightByNormalAndIn(RootedBinaryTreeMapping firstMapping,
			RootedBinaryTreeMapping inMapping, DTLReconciliation reconciliation) {
		BigInteger inCount = new BigInteger("0");
		
		for(Node inNode: inMapping.getInSpeciesNodes()) {
			Node geneNode = inMapping.geneNode;
			RootedBinaryTreeMapping sampleInMapping = reconciliation.getMappingFromNodes(geneNode, inNode);
		inCount = 	inCount.add(sampleInMapping.numberOfSolutions);
		}
		optionWeights.add(firstMapping.numberOfSolutions.multiply(inCount));
		
		
	}

	public void addWeightByNormalAndNormal(RootedBinaryTreeMapping firstMapping,
			RootedBinaryTreeMapping secondMapping ) {
		optionWeights.add(firstMapping.numberOfSolutions.multiply(secondMapping.numberOfSolutions));
		
	}

	public void sumUpNumberOfOptions() {
		BigInteger sum = new BigInteger("0");
		for(BigInteger weight: optionWeights) {
			sum = sum.add(weight);
		}
		numberOfSolutions = sum;
	}

	public RootedBinaryTreeMapping getInMappingGreedily(DTLReconciliation reconciliation) {

		if(inSpeciesNodes.size()==0)return null;
		

		RootedBinaryTreeMapping greedyMapping = null;
		BigInteger numberOfSolutions  = new BigInteger("0");
		
		for(Node inNode : inSpeciesNodes )
		{
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, inNode);
			int comparison = mapping.getNumberOfSolutions().compareTo(numberOfSolutions);
			if(comparison==1) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
			}
			else if (comparison==0) {
				int choice =reconciliation.binaryChoice(); 
				if(choice == 1) {
					numberOfSolutions = mapping.getNumberOfSolutions();
					greedyMapping = mapping;
				}
			}
			// here i am
		//	if(inNode.isCandidadateNode(geneNode)) {
			//	inNode.addCandidateNode(geneNode);
		//	inNode.incrementMappingCount();
		//	}
		}
			return greedyMapping;
		
	
	}

	public RootedBinaryTreeMapping getOutMappingGreedily(DTLReconciliation reconciliation) {

		if(outSpeciesNodes.size()==0)return null;
		

		RootedBinaryTreeMapping greedyMapping = null;
		BigInteger numberOfSolutions  = new BigInteger("0");
		
		for(Node outNode : outSpeciesNodes )
		{
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, outNode);
			int comparison = mapping.getNumberOfSolutions().compareTo(numberOfSolutions);
			if(comparison==1) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
			}
			else if (comparison==0) {
				int choice =reconciliation.binaryChoice(); 
				if(choice == 1) {
					numberOfSolutions = mapping.getNumberOfSolutions();
					greedyMapping = mapping;
				}
			}
	//		outNode.incrementMappingCount();
		}
			return greedyMapping;
		
	
	}

	public RootedBinaryTreeMapping getInMappingUniformRandomly(DTLReconciliation reconciliation) {


		if(inSpeciesNodes.size()==0)return null;
		RootedBinaryTreeMapping uniformMapping = null;
		BigInteger numberOfSolutions  = new BigInteger("0");
	
		
		for (Node inNode : inSpeciesNodes)
		{
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, inNode);
			numberOfSolutions =  numberOfSolutions.add(mapping.getNumberOfSolutions());
		}
		double random = Math.random() * numberOfSolutions.doubleValue();
		for (int i = 0; i < inSpeciesNodes.size(); ++i)
		{
			Node inNode = inSpeciesNodes.get(i);
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, inNode);
		    random -= mapping.getNumberOfSolutions().doubleValue();
		    if (random <= 0.0d)
		    {
		        uniformMapping = mapping;
		        break;
		    }
		}
		return uniformMapping;	
	}

	public RootedBinaryTreeMapping getOutMappingUniformRandomly(DTLReconciliation reconciliation) {


		if(outSpeciesNodes.size()==0)return null;
		RootedBinaryTreeMapping uniformMapping = null;
		BigInteger numberOfSolutions  = new BigInteger("0");
	
		
		for (Node inNode : outSpeciesNodes)
		{
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, inNode);
			numberOfSolutions =  numberOfSolutions.add(mapping.getNumberOfSolutions());
		}
		double random = Math.random() * numberOfSolutions.doubleValue();
		for (int i = 0; i < outSpeciesNodes.size(); ++i)
		{
			Node outNode = outSpeciesNodes.get(i);
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, outNode);
		    random -= mapping.getNumberOfSolutions().doubleValue();
		    if (random <= 0.0d)
		    {
		        uniformMapping = mapping;
		        break;
		    }
		}
		return uniformMapping;	
	}

	public boolean isTransferPossible() {
		if(cost== transferCost)return true;
		return false;
	}

	public RootedBinaryTreeMapping getTransferBiasedInMapping(DTLReconciliation reconciliation) {


		if(inSpeciesNodes.size()==0)return null;
		

		RootedBinaryTreeMapping greedyMapping = null;
		BigInteger numberOfSolutions  = new BigInteger("0");
		
		for(Node inNode : inSpeciesNodes )
		{
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, inNode);
			int comparison = mapping.getNumberOfSolutions().compareTo(numberOfSolutions);
			if(mapping.isTransferPossible()) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
				break;
			}
			else	if(comparison==1) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
			}
			else if (comparison==0) {
				int choice =reconciliation.binaryChoice(); 
				if(choice == 1) {
					numberOfSolutions = mapping.getNumberOfSolutions();
					greedyMapping = mapping;
				}
			}
			// here i am
				inNode.addCandidateNode(geneNode);
			inNode.incrementMappingCount();
			
		}
			return greedyMapping;
		
	
	
	}

	public RootedBinaryTreeMapping getTransferBiasedOutMapping(DTLReconciliation reconciliation) {


		if(outSpeciesNodes.size()==0)return null;
		

		RootedBinaryTreeMapping greedyMapping = null;
		BigInteger numberOfSolutions  = new BigInteger("0");
		
		for(Node outNode : outSpeciesNodes )
		{
			RootedBinaryTreeMapping mapping = reconciliation.getMappingFromNodes(geneNode, outNode);
			int comparison = mapping.getNumberOfSolutions().compareTo(numberOfSolutions);
			if(mapping.isTransferPossible()) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
				break;
			}
			else if(comparison==1) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
			}
			else if (comparison==0) {
				int choice =reconciliation.binaryChoice(); 
				if(choice == 1) {
					numberOfSolutions = mapping.getNumberOfSolutions();
					greedyMapping = mapping;
				}
			}
		}
			return greedyMapping;
		
	
	
	}

	public int getGreedyOption() {
		BigInteger greedyOptionValue  = new BigInteger("0");
		int greedyOption = 0;
		int i = 0;
		for(BigInteger weight : optionWeights )
		{

			int comparison = weight.compareTo(greedyOptionValue);
			
			 if(comparison==1) {
				greedyOptionValue = weight;
				greedyOption = options.get(i);
			}
			else if (comparison==0) {
				int choice = new Random().nextInt(2);
				if(choice == 1) {
					greedyOptionValue = weight;
					greedyOption = options.get(i);
				}
			}
			 i++;
		}
		return greedyOption;
	}

	public void updateMappingAndDuplicationCounts() {
		//if(speciesNode.isCandidadateNode(geneNode)) {
			speciesNode.addCandidateNode(geneNode);
			speciesNode.incrementMappingCount();
	//	}
		for(int option: options) {
			if(option== 1 || (option>5 && option<15)) {
				speciesNode.incrementDuplicationCount();
			}
		}
	}

	
}
