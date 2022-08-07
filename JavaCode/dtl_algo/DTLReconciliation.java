package dtl_algo;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.*;
import mapping.*;
import node.GeneNodeRooted;
import node.Node;
import tree.*;
import utility.FileUtility;

public class DTLReconciliation {
	SpeciesTreeRooted speciesTree; // species tree instance
	GeneTreeRooted geneTree; // gene tree instance
	public static int DUPLICATION = 2; // cost of duplication
	public static int LOSS = 1; // cost of loss
	public static int TRANSFER = 3; // cost of loss
	public static int ADDITIVE_TRANSFER = 4; // cost of additive transfer
	public static int REPLACEMENT_TRANSFER = 3; // cost of additive transfer
	int reconciliationCost; // final cost of reconciliation
	HashMap<Node, ArrayList<RootedBinaryTreeMapping>> optimalMappings = new HashMap<>(); // list of all optimal mappings
	RootedBinaryTreeMapping[][] allPossibleMappings;
	boolean inputError = false;
	int numberOfTransfer = 0;
	int numberOfDuplication = 0;
	// setters and getters of all species tree and gene tree
	public SpeciesTreeRooted getSpeciesTree() {
		return speciesTree;
	}

	public void setSpeciesTree(SpeciesTreeRooted speciesTree) {
		this.speciesTree = speciesTree;
	}

	public GeneTreeRooted getGeneTree() {
		return geneTree;
	}

	public void setGeneTree(GeneTreeRooted geneTree) {
		this.geneTree = geneTree;
	}

	public boolean isInputError() {
		return inputError;
	}

	public void setInputError(boolean inputError) {
		this.inputError = inputError;
	}

	// main function
	public static void main(String[] args) throws IOException {

		 processSingleInputOutput(args);
		//createMLInput(args);
		// createRealDataInput(args);
	}

	private static void createRealDataInput(String[] args) throws IOException {

		String inputFileName = "";
		String outputFileName = "";
		String trueLabelFileName = "";

		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-i")) {
				i++;
				inputFileName = args[i];
			} else if (args[i].equals("-o")) {
				i++;
				outputFileName = args[i];
			} else if (args[i].equals("-L")) {
				i++;
				LOSS = Integer.parseInt(args[i]);
			} else if (args[i].equals("-D")) {
				i++;
				DUPLICATION = Integer.parseInt(args[i]);
			} else if (args[i].equals("-T")) {
				i++;
				ADDITIVE_TRANSFER = Integer.parseInt(args[i]);
			} else if (args[i].equals("-R")) {
				i++;
				REPLACEMENT_TRANSFER = Integer.parseInt(args[i]);
			} else if (args[i].equals("-real")) {
				i++;
				trueLabelFileName = args[i];
			}
		}
	
		DTLReconciliation reconciliation = new DTLReconciliation();
		FileUtility.processInput(inputFileName, reconciliation);
		
		if (reconciliation.isInputError())
			return;
		FileWriter writer = new FileWriter(outputFileName, true);
		reconciliation.performPreprocessings();
		reconciliation.performAllDTLReconciliation();

		reconciliation.sampleGreedyOptionSolution();
		//System.out.println(inputFileName);
		reconciliation.printInputForRealData(writer,  trueLabelFileName);
		writer.close();

	}

	private void printInputForRealData(FileWriter writer, String trueLabelFileName) throws IOException {
		HashMap<String, Integer> currentHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> lossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> mappingCountHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> siblingHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> singleBranchLossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> tripleBranchLossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> anyBranchLossHeuristicMap = new HashMap<String, Integer>();

		updateHeuristicMaps(currentHeuristicMap, lossHeuristicMap, mappingCountHeuristicMap,
				siblingHeuristicMap, singleBranchLossHeuristicMap, tripleBranchLossHeuristicMap,
				anyBranchLossHeuristicMap);

		ArrayList<Node> allInternalGeneNodes = geneTree.getAllPostOrderInternalNodes();
		for (Node geneNode : allInternalGeneNodes) {
			String name = geneNode.getInternalNodeName();
			int trueLabel = 2;
			int currentHeuristicLabel = 2;
			int lossHeuristicLabel = 2;
			int mappingCountHeuristicLabel = 2;
			
			int siblingHeuristicLabel = 2;
			int singleLevelLossHeuristicLabel = 2;
			int tripleLevelLossHeuristicLabel = 2;
			int anyLevelLossHeuristicLabel = 2;
			int distanceOfTransfer = 0;
			if (!currentHeuristicMap.containsKey(name)  ) {
				continue;
			}
			else if(!isLeafToLeafTransfer(geneNode)){
				continue;
			}
			//within = 1 accross = 0=========================== =========== 
			 if(getWithinOrAccrossTransferValue(geneNode)==1)
				continue;
			else {
				distanceOfTransfer = getTransferDistance(geneNode );
				currentHeuristicLabel = currentHeuristicMap.get(name);
				lossHeuristicLabel = lossHeuristicMap.get(name);
				mappingCountHeuristicLabel = mappingCountHeuristicMap.get(name);
				siblingHeuristicLabel = siblingHeuristicMap.get(name);
				 singleLevelLossHeuristicLabel = singleBranchLossHeuristicMap.get(name);
				 
				 tripleLevelLossHeuristicLabel = tripleBranchLossHeuristicMap.get(name);
				 anyLevelLossHeuristicLabel = anyBranchLossHeuristicMap.get(name);
			}

			int height = geneNode.getHeight();
			//int depth = geneNode.getDepth();
		//	int numberOfDescendants = geneNode.getDescendents().size();
			//within = 1 accross = 0=========================== =========== 
		//	trueLabel = getWithinOrAccrossTransferValue(geneNode);
		//	if(geneNode.getName().isEmpty())geneNode.setName("m1");
			//original version//
			/*
			String data = geneNode.getName() + "," + height + "," + depth + "," + numberOfDescendants + "," + mappingCountHeuristicLabel + ","
					+ siblingHeuristicLabel  + "," + lossHeuristicLabel + "," + singleLevelLossHeuristicLabel 
					+","+tripleLevelLossHeuristicLabel+","+anyLevelLossHeuristicLabel+","
					+ currentHeuristicLabel +","+ distanceOfTransfer + ","+trueLabel + "\n";
			*/
			GeneNodeRooted rootedGeneNode = (GeneNodeRooted) geneNode;
			RootedBinaryTreeMapping mapping = rootedGeneNode.getMapping();
			Node transferMapping = mapping.getSpeciesNode();
			Node transferRecipient = mapping.getAdditiveTransferRecipient();
			//Lina kloub version
			String data = geneNode.getInternalNodeName() + "," + height  + "," + mappingCountHeuristicLabel + ","
					+ lossHeuristicLabel + "," + singleLevelLossHeuristicLabel 
					+","+tripleLevelLossHeuristicLabel+","
					+ currentHeuristicLabel +","+distanceOfTransfer+","+ trueLabel + "\n";
			writer.write(data);
		
		}
	}

	private int getWithinOrAccrossTransferValue(Node geneNode) {
		GeneNodeRooted rootedGeneNode = (GeneNodeRooted) geneNode;
		RootedBinaryTreeMapping mapping = rootedGeneNode.getMapping();
		Node transferMapping = mapping.getSpeciesNode();
		Node transferRecipient = mapping.getAdditiveTransferRecipient();
		String mappingName = transferMapping.getRealName();
		String recipientName = transferRecipient.getRealName();
		//System.out.println("Mapping: "+mappingName);
		//System.out.println("Recipient: "+ recipientName);
		if(!mappingName.contains("_") || !recipientName.contains("_"))
			return 0;
		else {
			String firstPartOfDonor = mappingName.split("_")[0];
			String firstPartOfRecipient = recipientName.split("_")[0];
			if(firstPartOfDonor.equals(firstPartOfRecipient))return 1;
			else return 0;
		}
		
	}

	private boolean isLeafToLeafTransfer(Node geneNode) {
		GeneNodeRooted rootedGeneNode = (GeneNodeRooted) geneNode;
		RootedBinaryTreeMapping mapping = rootedGeneNode.getMapping();
		Node transferMapping = mapping.getSpeciesNode();
		Node transferRecipient = mapping.getAdditiveTransferRecipient();
		if(transferRecipient.isLeaf() && transferMapping.isLeaf()) {
			return true;
		}
		return false;
	}

	private static void createMLInput(String[] args) throws IOException {
		String inputFileName = "";
		String outputFileName = "";
		String trueLabelFileName = "";
		// System.out.println(args);
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-i")) {
				i++;
				inputFileName = args[i];
			} else if (args[i].equals("-o")) {
				i++;
				outputFileName = args[i];
			} else if (args[i].equals("-L")) {
				i++;
				LOSS = Integer.parseInt(args[i]);
			} else if (args[i].equals("-D")) {
				i++;
				DUPLICATION = Integer.parseInt(args[i]);
			} else if (args[i].equals("-T")) {
				i++;
				ADDITIVE_TRANSFER = Integer.parseInt(args[i]);
			} else if (args[i].equals("-R")) {
				i++;
				REPLACEMENT_TRANSFER = Integer.parseInt(args[i]);
			} else if (args[i].equals("-real")) {
				i++;
				trueLabelFileName = args[i];
			}
		}
		DTLReconciliation reconciliation = new DTLReconciliation();
		FileUtility.processInput(inputFileName, reconciliation);
		if (reconciliation.isInputError())
			return;
		FileWriter writer = new FileWriter(outputFileName, true);
		reconciliation.performPreprocessings();
		reconciliation.performAllDTLReconciliation();

		reconciliation.sampleGreedyOptionSolution();
		// reconciliation.sampleUniformRandomSolution();
		// reconciliation.modifySolutionWithHeuristic();
		HashMap<String, Integer> trueLabelMap = getTrueLabelMap(trueLabelFileName);
		// reconciliation.printOutputOfDTRL(writer);

		reconciliation.printInputForSupervisedML(writer, trueLabelMap);
	//	reconciliation.getNumberOfTransfersAndDuplications();
	//	int numDuplication = reconciliation.numberOfDuplication;
	//	int numTransfers = reconciliation.numberOfTransfer;
		
		//FileWriter geneWriter = new FileWriter("geneTreeSize.txt",true);
		//geneWriter.write(reconciliation.geneTree.getLeafSet().size()+ " "+ numDuplication+" "+numTransfers+"\n");
	//	geneWriter.close();
		writer.close();
	}

	

	private void printInputForSupervisedML(FileWriter writer, HashMap<String, Integer> trueLabelMap)
			throws IOException {
		
		HashMap<String, Integer> currentHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> lossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> mappingCountHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> siblingHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> singleBranchLossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> tripleBranchLossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> anyBranchLossHeuristicMap = new HashMap<String, Integer>();

		updateHeuristicMaps(currentHeuristicMap, lossHeuristicMap, mappingCountHeuristicMap,
				siblingHeuristicMap, singleBranchLossHeuristicMap, tripleBranchLossHeuristicMap,
				anyBranchLossHeuristicMap);

		ArrayList<Node> allInternalGeneNodes = geneTree.getAllPostOrderInternalNodes();
		//System.out.println("Debugging inside all nodes");
		for (Node geneNode : allInternalGeneNodes) {
			String name = geneNode.getInternalNodeName();
			int trueLabel = 2;
			int currentHeuristicLabel = 2;
			int lossHeuristicLabel = 2;
			int mappingCountHeuristicLabel = 2;
		
			int siblingHeuristicLabel = 2;
			int singleLevelLossHeuristicLabel = 2;
			int tripleLevelLossHeuristicLabel = 2;
			int anyLevelLossHeuristicLabel = 2;
			
			if (!trueLabelMap.containsKey(name) /*&& !currentHeuristicMap.containsKey(name)*/ ) {
				continue;
			}
			if (trueLabelMap.containsKey(name)) {
				trueLabel = trueLabelMap.get(name);

			}
			
			if (currentHeuristicMap.containsKey(name)) {	
				currentHeuristicLabel = currentHeuristicMap.get(name);
				lossHeuristicLabel = lossHeuristicMap.get(name);
				mappingCountHeuristicLabel = mappingCountHeuristicMap.get(name);
				siblingHeuristicLabel = siblingHeuristicMap.get(name);
				singleLevelLossHeuristicLabel = singleBranchLossHeuristicMap.get(name);
				 
				 tripleLevelLossHeuristicLabel = tripleBranchLossHeuristicMap.get(name);
				 anyLevelLossHeuristicLabel = anyBranchLossHeuristicMap.get(name);
			}
			else continue;

			int height = geneNode.getHeight();
			
			String data = geneNode.getInternalNodeName() + "," + height  + "," + mappingCountHeuristicLabel + ","
					+ lossHeuristicLabel + "," + singleLevelLossHeuristicLabel 
					+","+tripleLevelLossHeuristicLabel+","
					+ currentHeuristicLabel +","+ trueLabel + "\n";

			writer.write(data);
	
		}
	}
	private int getTransferDistance(Node geneNode) {
		GeneNodeRooted rootedGeneNode = (GeneNodeRooted) geneNode;
		
		RootedBinaryTreeMapping mapping = rootedGeneNode.getMapping();
		
		Node mappedNode = mapping.getSpeciesNode();
		
		Node transferRecipient = mapping.getAdditiveTransferRecipient();
		return mappedNode.getDistance(transferRecipient);
	}





	
//Abhijit
	private void updateHeuristicMaps(HashMap<String, Integer> currentHeuristicMap,
			HashMap<String, Integer> lossHeuristicMap, HashMap<String, Integer> mappingCountHeuristicMap,
			  HashMap<String, Integer> siblingHeuristicMap,
			HashMap<String, Integer> singleBranchLossHeuristicMap,
			HashMap<String, Integer> tripleBranchLossHeuristicMap, HashMap<String, Integer> anyBranchLossHeuristicMap) {

		ArrayList<Node> allPreOrderInternalGeneNodes = geneTree.getAllPreOrderInternalNodes();
		ArrayList<Node> geneLeaves = geneTree.getLeafSet();
		ArrayList<Node> speciesLeaves = speciesTree.getLeafSet();
		HashMap<String, Integer> actualGeneFrequencies = getFrequencyMapFromNodes(geneLeaves);
		HashMap<String, Integer> inferredGeneFrequencies = new HashMap<String, Integer>();

		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node speciesNode = mapping.getSpeciesNode();

			if (geneNode.getEvent().equals("Duplication")) {
				speciesNode.incrementDuplicationCount();
			} else if (geneNode.getEvent().equals("Transfer")) {

				speciesNode.incrementTransferCount();
				Node transferRecipient = mapping.getAdditiveTransferRecipient();
				transferRecipient.incrementTransferRecipientCount();
			}
		}
		initializeInferredFrequencies(inferredGeneFrequencies, speciesLeaves);
		updateMappingAndDuplicationCounts();
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();

			if (geneNode.getEvent().equals("Transfer")) {
				Node childCausingTransfer = null;
				if(geneNode.getInternalNodeName()==null) {
					geneNode.setInternalNodeName(geneNode.getName());
				}
				if (geneNode.getOption() == 2 || geneNode.getOption() == 15) {
					childCausingTransfer = geneNode.getRight();
				} else {
					childCausingTransfer = geneNode.getLeft();
				}
				Node additiveRecipient = mapping.getAdditiveTransferRecipient();

				if (additiveRecipient.additiveByMappingCount(geneNode)) {
					mappingCountHeuristicMap.put(geneNode.getInternalNodeName(), 1);
				//	System.out.println("found");
				} else {

					mappingCountHeuristicMap.put(geneNode.getInternalNodeName(), 0);
				}
				if (!additiveRecipient.anyAncestorBranchHavingLosses()) {
					lossHeuristicMap.put(geneNode.getInternalNodeName(), 1);
				} else {
					lossHeuristicMap.put(geneNode.getInternalNodeName(), 0);

				}
				if (!currentHeuristic(additiveRecipient, mapping, actualGeneFrequencies, inferredGeneFrequencies)) {
					currentHeuristicMap.put(geneNode.getInternalNodeName(), 1);
				} else {
					currentHeuristicMap.put(geneNode.getInternalNodeName(), 0);

				}
				

				
				if (!additiveRecipient.anyAncestorBranchHavingLosses(1)) {
					singleBranchLossHeuristicMap.put(geneNode.getInternalNodeName(), 1);
				} else {
					singleBranchLossHeuristicMap.put(geneNode.getInternalNodeName(), 0);

				}
				if (!additiveRecipient.anyAncestorBranchHavingLosses(3)) {
					tripleBranchLossHeuristicMap.put(geneNode.getInternalNodeName(), 1);
				} else {
					tripleBranchLossHeuristicMap.put(geneNode.getInternalNodeName(), 0);

				}
				if (!additiveRecipient.anyAncestorBranchHavingLosses(-1)) {
					anyBranchLossHeuristicMap.put(geneNode.getInternalNodeName(), 1);
				} else {
					anyBranchLossHeuristicMap.put(geneNode.getInternalNodeName(), 0);

				}
			}
		}

	}

	
	private static HashMap<String, Integer> getTrueLabelMap(String trueLabelFileName) throws FileNotFoundException {
		// System.out.println(trueLabelFileName);
		Scanner scanner = new Scanner(new File(trueLabelFileName));
		HashMap<String, Integer> map = new HashMap<>();

		for (int i = 0; i < 3; i++)
			scanner.nextLine();
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			line = line.replaceAll("\t", " ");
			line = line.replaceAll("( )+", " ");
			// System.out.println(line);
			String[] lineParts = line.split(" ");
			String[] nodeNameParts = lineParts[0].split("_");
			String nodeName = "H" + nodeNameParts[1] + "_" + nodeNameParts[0].substring(1);
			String event = lineParts[2];
			if (event.contains("ADDITIVE_TRANSFER"))
				map.put(nodeName, 1);
			else if (event.contains("REPLACING_TRANSFER"))
				map.put(nodeName, 0);
			
		}
		scanner.close();
		return map;
	}

	// takes input folder name and output folder name and does DTL reconciliation
	// with heuristic, outputs only cost though
	public static void processFolderInputOutputOnlyCost(String inputFolderName, String outputFolderName)
			throws IOException {
		File folder = new File(inputFolderName);
		File[] listOfFiles = folder.listFiles();

		for (int i = 0; i < listOfFiles.length; i++) {
			DTLReconciliation reconciliation = new DTLReconciliation();
			if (listOfFiles[i].isFile()) {

				String inputFileName = listOfFiles[i].getName();
				FileUtility.processInput(inputFolderName + "\\" + inputFileName, reconciliation);
				if (reconciliation.isInputError())
					continue;
				reconciliation.performPreprocessings();

				reconciliation.performDTLReconciliationWithOptionTracking();
				reconciliation.sampleRandomSolutionWithLossTracking();

				reconciliation.performDTLReplacementHeuristic();
				PrintWriter writer = new PrintWriter(new File(outputFolderName + inputFileName));

				writer.println(reconciliation.reconciliationCost);
				writer.println(reconciliation.geneTree.getLeafSet().size());
				writer.close();
			}
		}
	}

	// takes single input file name and single output filename and does basic DTL
	// reconciliation
	public static void processSingleInputOutput(String[] args) throws IOException {
		String inputFileName = "";
		String outputFileName = "";
		String heuristic = "ml";
		for (int i = 0; i < args.length; i++) {
			if (args[i].equals("-i")) {
				i++;
				inputFileName = args[i];
			} else if (args[i].equals("-o")) {
				i++;
				outputFileName = args[i];
			} else if (args[i].equals("-L")) {
				i++;
				LOSS = Integer.parseInt(args[i]);
			} else if (args[i].equals("-D")) {
				i++;
				DUPLICATION = Integer.parseInt(args[i]);
			} else if (args[i].equals("-T")) {
				i++;
				ADDITIVE_TRANSFER = Integer.parseInt(args[i]);
			} else if (args[i].equals("-R")) {
				i++;
				REPLACEMENT_TRANSFER = Integer.parseInt(args[i]);
			}
			else if (args[i].equals("-h")) {
				i++;
				heuristic = args[i];
			}
		}
		DTLReconciliation reconciliation = new DTLReconciliation();
		FileUtility.processInput(inputFileName, reconciliation);
		if (reconciliation.isInputError())
			return;
		
		reconciliation.performPreprocessings();
		reconciliation.performAllDTLReconciliation();

		reconciliation.sampleGreedyOptionSolution();
		//reconciliation.sampleUniformRandomSolution();
		if(heuristic.equals("geneFrequency")) {
			reconciliation.modifySolutionWithLeafCountHeuristic();
		}
		else if(heuristic.equals("mappingCount")) {
			reconciliation.modifySolutionWithMappingCountHeuristic();
		}
		else if(heuristic.equals("ml")) {
			
			reconciliation.modifySolutionByML(inputFileName);
		}
		else if(heuristic.equals("lostGene2")) {
		
			reconciliation.modifySolutionWithLossHeuristic();
		}
		else if(heuristic.equals("lostGene1")) {
			
			reconciliation.modifySolutionWithLossHeuristic(1);
		}
		else if(heuristic.equals("lostGene3")) {
			
			reconciliation.modifySolutionWithLossHeuristic(3);
		}
		else {
			System.out.println("Please provide a valid heuristic name for transfer classification: \n"
					+ "\t 'ml': Machine Learning method"			
					+ "\n\t 'lostGene1': Lost-gene(h=1) heuristic"
					+ "\n\t 'lostGene2': Lost-gene(h=2) heuristic"
					+ "\n\t 'lostGene3': Lost-gene(h=3) heuristic"
					+ "\n\t 'mappingCount': Mapping-count heuristic"
					+ "\n\t 'geneFrequency': Gene-frequency heuristic");
		}
		if(!outputFileName.isEmpty())
		{
			PrintWriter writer = new PrintWriter(new File(outputFileName));
			reconciliation.printOutputOfDTRL(writer);
			writer.close();
		}
		else
			reconciliation.printOutputOfDTRLConsole();
	
	}

	private void modifySolutionWithLossHeuristic(int i) {
		ArrayList<Node> allPreOrderInternalGeneNodes = geneTree.getAllPreOrderInternalNodes();
		ArrayList<Node> speciesLeaves = speciesTree.getLeafSet();
		HashMap<String, Integer> inferredGeneFrequencies = new HashMap<String, Integer>();
		
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node speciesNode = mapping.getSpeciesNode();

			if (geneNode.getEvent().equals("Duplication")) {
			//	speciesNode.incrementDuplicationCount();
			} else if (geneNode.getEvent().equals("Transfer")) {

				speciesNode.incrementTransferCount();
				Node transferRecipient = mapping.getAdditiveTransferRecipient();
				transferRecipient.incrementTransferRecipientCount();
			}
		}
		updateMappingAndDuplicationCounts();
		initializeInferredFrequencies(inferredGeneFrequencies, speciesLeaves);
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			
			if (geneNode.getEvent().equals("Transfer")) {
				Node additiveRecipient = mapping.getAdditiveTransferRecipient();
				
				if (additiveRecipient.anyAncestorBranchHavingLosses(i)) {
					geneNode.setEvent("Replacing Transfer");
				} else {
					geneNode.setEvent("Additive Transfer");

				}
			}
		}
	}

	private void printOutputOfDTRLConsole() {
		//System.out.println("Trying to print");
				System.out.println("\n\n ------------ Reconciliation for Gene Tree 1 (rooted) -------------");
				System.out.println("Species Tree:");
				System.out.println(speciesTree.getNewickOutputFormat() + "\n");
				System.out.println("Gene Tree:");
				System.out.println(geneTree.getNewickOutputFormat() + "\n");
				System.out.println("Reconciliation:");
				ArrayList<Node> allPostOrderGeneNodes = geneTree.getAllNodesByPostOrder();
				int numberOfDuplication = 0;
				int numberOfTransfer = 0;

				int numberOfLoss = 0;
				for (Node node : allPostOrderGeneNodes) {
					GeneNodeRooted geneNode = (GeneNodeRooted) node;
					if (geneNode.isLeaf()) {
						System.out.println(geneNode.getFullName() + ": " + "Leaf Node");
						continue;
					}
					// writer.println("Internal node name: "+node.getInternalNodeName());
					System.out.print(geneNode.getName() + " = " + geneNode.getLCAFormat() + ": " + geneNode.getEvent()
							+ ", Mapping --> " + geneNode.getMappedSpeciesNode().getName());

					if (geneNode.getEvent().contains("Transfer")) {
						System.out.print(", Recipient --> " + geneNode.getMapping().getAdditiveTransferRecipient());
						numberOfTransfer++;

					}
					if (geneNode.getEvent().equals("Duplication"))
						numberOfDuplication++;

					if (geneNode.getMapping().getLosses().size() > 0) {
						numberOfLoss += geneNode.getMapping().getLostEdges().size();
					}
					System.out.println();

				}
				numberOfLoss = reconciliationCost - numberOfDuplication*DUPLICATION - numberOfTransfer*TRANSFER;
				System.out.print("\nThe minimum reconciliation cost is: " + reconciliationCost);
				System.out.println(" (Duplications: " + numberOfDuplication + ", Transfers: " + numberOfTransfer + ", Losses: "
						+ numberOfLoss + ")");
				System.out.println("Total number of optimal solutions: " + getNumberOfSolutions());
				System.out.println("Total number of candidates for gene birth: " + optimalMappings.get(geneTree.getRoot()).size());
				System.out.close();
	}

	private void modifySolutionByML(String inputFileName) throws FileNotFoundException {
		
		
		HashMap<String, Integer> currentHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> lossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> mappingCountHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> siblingHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> singleBranchLossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> tripleBranchLossHeuristicMap = new HashMap<String, Integer>();
		HashMap<String, Integer> anyBranchLossHeuristicMap = new HashMap<String, Integer>();

		updateHeuristicMaps(currentHeuristicMap, lossHeuristicMap, mappingCountHeuristicMap,
				siblingHeuristicMap, singleBranchLossHeuristicMap, tripleBranchLossHeuristicMap,
				anyBranchLossHeuristicMap);

		ArrayList<Node> allInternalGeneNodes = geneTree.getAllPostOrderInternalNodes();
		// get a random string of length 7
		String randomString = getRanDomString(7);
		String pythonInputFileName = inputFileName+"_"+randomString+".txt";
		PrintWriter pythonInputWriter = new PrintWriter(new File(pythonInputFileName)); 
		pythonInputWriter.write("id,height,mappingCountHeuristic,lostGene2Heuristic,lostGene1Heuristic,lostGene3Heuristic,"
				+ "geneFrequencyHeuristic,"  + "trueLabel\r\n");
		int transferCount = 0;
		for (Node node : allInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			String name = geneNode.getInternalNodeName();
			int trueLabel = 2;
			int currentHeuristicLabel = 2;
			int lossHeuristicLabel = 2;
			int mappingCountHeuristicLabel = 2;
		
			
			int singleLevelLossHeuristicLabel = 2;
			int tripleLevelLossHeuristicLabel = 2;
			
			if (currentHeuristicMap.containsKey(name)) {	
				currentHeuristicLabel = currentHeuristicMap.get(name);
				lossHeuristicLabel = lossHeuristicMap.get(name);
				mappingCountHeuristicLabel = mappingCountHeuristicMap.get(name);
			
				singleLevelLossHeuristicLabel = singleBranchLossHeuristicMap.get(name);
				 
				 tripleLevelLossHeuristicLabel = tripleBranchLossHeuristicMap.get(name);
				
			}
			else continue;

			int height = geneNode.getHeight();
			
			if (geneNode.getEvent().equals("Transfer")) {
				transferCount++;
				String data = geneNode.getInternalNodeName() + "," + height + "," +     mappingCountHeuristicLabel + ","
						+  lossHeuristicLabel + "," + singleLevelLossHeuristicLabel 
						+","+tripleLevelLossHeuristicLabel+","
						+ currentHeuristicLabel +","+ trueLabel + "\n";
				pythonInputWriter.write(data);
			}
		}
		
		pythonInputWriter.close();
		ArrayList<Integer> mlOutputs = new ArrayList<Integer>();
		 String outputOfPython = null;

	        try {
	            Process p = Runtime.getRuntime().exec("python DTRL_Classifier.py "+pythonInputFileName);
	            
	            BufferedReader stdInput = new BufferedReader(new 
	                 InputStreamReader(p.getInputStream()));

	            BufferedReader stdError = new BufferedReader(new 
	                 InputStreamReader(p.getErrorStream()));

	            // read the output from the command
	           // System.out.println("Here is the standard output of the command:\n");
	            while ((outputOfPython = stdInput.readLine()) != null) {
	            
	            	//System.out.println(outputOfPython);
	                
	                
	                for(int i = 0 ;i<outputOfPython.length();i++)
	                {
	                	char c = outputOfPython.charAt(i);
	                	if(c=='0')
	                		mlOutputs.add(0);
	                	else
	                		mlOutputs.add(1);
	                		
	                }
	            }
	       
	            //System.out.println("Here is the standard error of the command (if any):\n");
	            while ((outputOfPython = stdError.readLine()) != null) {
	                System.out.println(outputOfPython);
	            }
	           
	            int indexOfOutput = 0;
	            for (Node node : allInternalGeneNodes) {
	            	GeneNodeRooted geneNode = (GeneNodeRooted) node;
	            	
	            	if (geneNode.getEvent().equals("Transfer")) {
	            		if(mlOutputs.get(indexOfOutput)==1) {
	            			geneNode.setEvent("Additive Transfer");
	            		}
	            		else 
	            			geneNode.setEvent("Replacing Transfer");
	            		indexOfOutput++;
	            	}
	            	
	            }
	          
	        }
	        catch (IOException e) {
	            //System.out.println("exception happened - here's what I know: ");
	            e.printStackTrace();
	            System.exit(-1);
	        }
	        
	}

	private String getRanDomString(int length) {
		String candidateChars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
		StringBuilder sb = new StringBuilder();
	    Random random = new Random();
	    for (int i = 0; i < length; i++) {
	        sb.append(candidateChars.charAt(random.nextInt(candidateChars
	                .length())));
	    }

	    return sb.toString();
		
	}

	private void sampleUniformRandomSolution() {

		ArrayList<Node> allPreOrderGeneNodes = geneTree.getAllNodesByPreOrder();
		Random random = new Random();
		for (Node geneNode : allPreOrderGeneNodes) {
			RootedBinaryTreeMapping optimalMapping = null;
			ArrayList<RootedBinaryTreeMapping> optimalMappingsOfGeneNode = optimalMappings.get(geneNode);
			GeneNodeRooted geneNodeRooted = (GeneNodeRooted) geneNode;
			if (geneNode.isRoot()) {

				optimalMapping = getUniformRandomSolution(optimalMappingsOfGeneNode);// getGreedyOptionSolution(optimalMappingsOfGeneNode);

			} else {
				GeneNodeRooted parentNode = (GeneNodeRooted) geneNodeRooted.getParent();

				RootedBinaryTreeMapping parentMapping = parentNode.getMapping();
				Node parentImage = parentMapping.getSpeciesNode();
				Node rightOfParentImage = parentMapping.getSpeciesNode().getRight();
				Node leftOfParentImage = parentMapping.getSpeciesNode().getLeft();

				int parentOption = parentNode.getOption();
				if (parentOption == 1) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 2 || parentOption == 15) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, parentImage);
						parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());
					} else {
						optimalMapping = getUniformOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 3 || parentOption == 16) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, parentImage);
						parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());

					} else {
						optimalMapping = getUniformOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 4) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, leftOfParentImage);

					} else {

						optimalMapping = getUniformRandomInMappingFromPair(geneNode, rightOfParentImage);

					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 5) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, rightOfParentImage);
					} else {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, leftOfParentImage);
					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 6) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());
					} else {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, rightOfParentImage);// inSearchWithLossTracking(optimalMappingsOfGeneNode,
						// rightOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 7) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 8) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());

					} else {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 9) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 10) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 11) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());
					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 12) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getUniformRandomInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());

					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 13) {
					optimalMapping = getUniformRandomInMappingFromPair(geneNode, leftOfParentImage);
					parentMapping.addLoss(parentImage.getRight());
					parentMapping.addLostBranch(parentImage.getRightBranch());
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
				} else if (parentOption == 14) {
					optimalMapping = getUniformRandomInMappingFromPair(geneNode, rightOfParentImage);
					parentMapping.addLoss(parentImage.getLeft());
					parentMapping.addLostBranch(parentImage.getLeftBranch());
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				}
			}

			geneNodeRooted.setMappedSpeciesNode(optimalMapping.getSpeciesNode());
			geneNodeRooted
					.setOption(optimalMapping.getOptions().get(random.nextInt(optimalMapping.getOptions().size())));
			geneNodeRooted.setMapping(optimalMapping);
			geneNodeRooted.setEventByOption();

		}

	}

	private RootedBinaryTreeMapping getUniformRandomInMappingFromPair(Node geneNode, Node speciesNode) {
		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);

		return mapping.getInMappingUniformRandomly(this);

	}

	private RootedBinaryTreeMapping getUniformRandomSolution(
			ArrayList<RootedBinaryTreeMapping> optimalMappingsOfGeneNode) {

		RootedBinaryTreeMapping uniformMapping = null;

		int minHeight = Integer.MAX_VALUE;
		for (RootedBinaryTreeMapping mapping : optimalMappingsOfGeneNode) {
			int height = mapping.getSpeciesNode().getHeight();
			//System.out.println(mapping.getSpeciesNode() + " " + height);
			if (height < minHeight) {
				uniformMapping = mapping;
				minHeight = height;
			}
		}
		return uniformMapping;

	}
	public void sampleTransferBiasedSolution() {


		ArrayList<Node> allPreOrderGeneNodes = geneTree.getAllNodesByPreOrder();
		for (Node geneNode : allPreOrderGeneNodes) {
			RootedBinaryTreeMapping optimalMapping = null;
			ArrayList<RootedBinaryTreeMapping> optimalMappingsOfGeneNode = optimalMappings.get(geneNode);
			GeneNodeRooted geneNodeRooted = (GeneNodeRooted) geneNode;
			if (geneNode.isRoot()) {

				optimalMapping = getTransferBiasedSolution(optimalMappingsOfGeneNode);

			} else {
				GeneNodeRooted parentNode = (GeneNodeRooted) geneNodeRooted.getParent();

				RootedBinaryTreeMapping parentMapping = parentNode.getMapping();
				Node parentImage = parentMapping.getSpeciesNode();
				Node rightOfParentImage = parentMapping.getSpeciesNode().getRight();
				Node leftOfParentImage = parentMapping.getSpeciesNode().getLeft();

				int parentOption = parentNode.getOption();
				if (parentOption == 1) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 2 || parentOption == 15) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, parentImage);
						parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());
					} else {
						optimalMapping = getTransferBiasedOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 3 || parentOption == 16) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, parentImage);
						parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());

					} else {
						optimalMapping = getTransferBiasedOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 4) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, leftOfParentImage);

					} else {

						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, rightOfParentImage);

					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 5) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, rightOfParentImage);
					} else {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, leftOfParentImage);
					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 6) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());
					} else {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, rightOfParentImage);// inSearchWithLossTracking(optimalMappingsOfGeneNode,
						// rightOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 7) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 8) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());

					} else {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 9) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 10) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 11) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());
					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 12) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getTransferBiasedInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());

					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 13) {
					optimalMapping = getTransferBiasedInMappingFromPair(geneNode, leftOfParentImage);
					parentMapping.addLoss(parentImage.getRight());
					parentMapping.addLostBranch(parentImage.getRightBranch());
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
				} else if (parentOption == 14) {
					optimalMapping = getTransferBiasedInMappingFromPair(geneNode, rightOfParentImage);
					parentMapping.addLoss(parentImage.getLeft());
					parentMapping.addLostBranch(parentImage.getLeftBranch());
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				}
			}

			geneNodeRooted.setMappedSpeciesNode(optimalMapping.getSpeciesNode());
			setTransferBiasedOption(geneNodeRooted, optimalMapping);
		
			geneNodeRooted.setMapping(optimalMapping);
			geneNodeRooted.setEventByOption();

		}

	
	}
	private void setTransferBiasedOption(GeneNodeRooted geneNodeRooted, RootedBinaryTreeMapping optimalMapping) {
		if(optimalMapping.getOptions().contains(2)) {
			geneNodeRooted.setOption(2);
		}
		else if(optimalMapping.getOptions().contains(3)) {
			geneNodeRooted.setOption(3);
		}
		else if(optimalMapping.getOptions().contains(15)) {
			geneNodeRooted.setOption(15);
		}
		else if(optimalMapping.getOptions().contains(16)) {
			geneNodeRooted.setOption(16);
		}
		else {
			Random random = new Random();
			geneNodeRooted.setOption(optimalMapping.getOptions().get(random.nextInt(optimalMapping.getOptions().size())));
		}
	}

	private RootedBinaryTreeMapping getTransferBiasedOutMappingFromPair(Node geneNode, Node speciesNode) {
		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);

		return mapping.getTransferBiasedOutMapping(this);
	}

	private RootedBinaryTreeMapping getTransferBiasedInMappingFromPair(Node geneNode, Node speciesNode) {
		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);

		return mapping.getTransferBiasedInMapping(this);
	}

	private void sampleGreedyOptionSolution() {

		ArrayList<Node> allPreOrderGeneNodes = geneTree.getAllNodesByPreOrder();
	
		for (Node geneNode : allPreOrderGeneNodes) {
			RootedBinaryTreeMapping optimalMapping = null;
			ArrayList<RootedBinaryTreeMapping> optimalMappingsOfGeneNode = optimalMappings.get(geneNode);
			GeneNodeRooted geneNodeRooted = (GeneNodeRooted) geneNode;
			if (geneNode.isRoot()) {

				optimalMapping = getGreedyOptionSolution(optimalMappingsOfGeneNode);

			} else {
				GeneNodeRooted parentNode = (GeneNodeRooted) geneNodeRooted.getParent();

				RootedBinaryTreeMapping parentMapping = parentNode.getMapping();
				Node parentImage = parentMapping.getSpeciesNode();
				Node rightOfParentImage = parentMapping.getSpeciesNode().getRight();
				Node leftOfParentImage = parentMapping.getSpeciesNode().getLeft();

				int parentOption = parentNode.getOption();
				if (parentOption == 1) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 2 || parentOption == 15) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getGreedyInMappingFromPair(geneNode, parentImage);
						parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());
					} else {
						optimalMapping = getGreedyOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());
					}
					
				} else if (parentOption == 3 || parentOption == 16) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getGreedyInMappingFromPair(geneNode, parentImage);
						parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());

					} else {
						optimalMapping = getGreedyOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());

					}
					
				} else if (parentOption == 4) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getGreedyInMappingFromPair(geneNode, leftOfParentImage);

					} else {

						optimalMapping = getGreedyInMappingFromPair(geneNode, rightOfParentImage);

					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 5) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getGreedyInMappingFromPair(geneNode, rightOfParentImage);
					} else {
						optimalMapping = getGreedyInMappingFromPair(geneNode, leftOfParentImage);
					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 6) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());
					} else {
						optimalMapping = getGreedyInMappingFromPair(geneNode, rightOfParentImage);// inSearchWithLossTracking(optimalMappingsOfGeneNode,
						// rightOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 7) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getGreedyInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 8) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());

					} else {
						optimalMapping = getGreedyInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 9) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getGreedyInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 10) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 11) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getGreedyInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getGreedyInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());
					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 12) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getGreedyInMappingFromPair(geneNode, rightOfParentImage);
						parentMapping.addLoss(parentImage.getRight());
						parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getGreedyInMappingFromPair(geneNode, leftOfParentImage);
						parentMapping.addLoss(parentImage.getLeft());
						parentMapping.addLostBranch(parentImage.getLeftBranch());

					}
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 13) {
					optimalMapping = getGreedyInMappingFromPair(geneNode, leftOfParentImage);
					parentMapping.addLoss(parentImage.getRight());
					parentMapping.addLostBranch(parentImage.getRightBranch());
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
				} else if (parentOption == 14) {
					optimalMapping = getGreedyInMappingFromPair(geneNode, rightOfParentImage);
					parentMapping.addLoss(parentImage.getLeft());
					parentMapping.addLostBranch(parentImage.getLeftBranch());
					parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
				}
			}

			geneNodeRooted.setMappedSpeciesNode(optimalMapping.getSpeciesNode());
			setOptionGreedily(geneNodeRooted, optimalMapping);
			geneNodeRooted.setMapping(optimalMapping);
			geneNodeRooted.setEventByOption();

		}

	}

	private void setOptionGreedily(GeneNodeRooted geneNodeRooted, RootedBinaryTreeMapping optimalMapping) {
		int greedyOption = optimalMapping.getGreedyOption();
		geneNodeRooted.setOption(greedyOption);
	}

	private RootedBinaryTreeMapping getUniformOutMappingFromPair(Node geneNode, Node speciesNode) {
		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);

		return mapping.getOutMappingUniformRandomly(this);
	}

	private RootedBinaryTreeMapping getGreedyOutMappingFromPair(Node geneNode, Node speciesNode) {
		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);

		return mapping.getOutMappingGreedily(this);
	}

	private RootedBinaryTreeMapping getGreedyInMappingFromPair(Node geneNode, Node speciesNode) {

		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);

		return mapping.getInMappingGreedily(this);

	}

	public int binaryChoice() {
		Random random = new Random();
		return random.nextInt(2);
	}

	private RootedBinaryTreeMapping getGreedyOptionSolution(
			ArrayList<RootedBinaryTreeMapping> optimalMappingsOfGeneNode) {
		RootedBinaryTreeMapping greedyMapping = null;
		BigInteger numberOfSolutions = new BigInteger("0");

		for (RootedBinaryTreeMapping mapping : optimalMappingsOfGeneNode) {
			int comparison = mapping.getNumberOfSolutions().compareTo(numberOfSolutions);
			if (comparison == 1) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
			} else if (comparison == 0) {
				int choice = binaryChoice();
				if (choice == 1) {
					numberOfSolutions = mapping.getNumberOfSolutions();
					greedyMapping = mapping;
				}
			}
		}
		return greedyMapping;
	}
	private RootedBinaryTreeMapping getTransferBiasedSolution(
			ArrayList<RootedBinaryTreeMapping> optimalMappingsOfGeneNode) {
		RootedBinaryTreeMapping greedyMapping = null;
		BigInteger numberOfSolutions = new BigInteger("0");

		for (RootedBinaryTreeMapping mapping : optimalMappingsOfGeneNode) {
			int comparison = mapping.getNumberOfSolutions().compareTo(numberOfSolutions);
			if(mapping.isTransferPossible()) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
				break;
			}
			else if (comparison == 1) {
				numberOfSolutions = mapping.getNumberOfSolutions();
				greedyMapping = mapping;
			} else if (comparison == 0) {
				int choice = binaryChoice();
				if (choice == 1 ) {
					numberOfSolutions = mapping.getNumberOfSolutions();
					greedyMapping = mapping;
				}
			}
		}
		return greedyMapping;
	}

	private void printOutputOfDTRL(PrintWriter writer) {
		//System.out.println("Trying to print");
		writer.println("\n\n ------------ Reconciliation for Gene Tree 1 (rooted) -------------");
		writer.println("Species Tree:");
		writer.println(speciesTree.getNewickOutputFormat() + "\n");
		writer.println("Gene Tree:");
		writer.println(geneTree.getNewickOutputFormat() + "\n");
		writer.println("Reconciliation:");
		ArrayList<Node> allPostOrderGeneNodes = geneTree.getAllNodesByPostOrder();
		int numberOfDuplication = 0;
		int numberOfTransfer = 0;

		int numberOfLoss = 0;
		for (Node node : allPostOrderGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			if (geneNode.isLeaf()) {
				writer.println(geneNode.getFullName() + ": " + "Leaf Node");
				continue;
			}
			// writer.println("Internal node name: "+node.getInternalNodeName());
			writer.print(geneNode.getName() + " = " + geneNode.getLCAFormat() + ": " + geneNode.getEvent()
					+ ", Mapping --> " + geneNode.getMappedSpeciesNode().getName());

			if (geneNode.getEvent().contains("Transfer")) {
				writer.print(", Recipient --> " + geneNode.getMapping().getAdditiveTransferRecipient());
				numberOfTransfer++;

			}
			if (geneNode.getEvent().equals("Duplication"))
				numberOfDuplication++;

			if (geneNode.getMapping().getLosses().size() > 0) {
				numberOfLoss += geneNode.getMapping().getLostEdges().size();
			}
			writer.println();

		}
		numberOfLoss = reconciliationCost - numberOfDuplication*DUPLICATION - numberOfTransfer*TRANSFER;
		writer.print("\nThe minimum reconciliation cost is: " + reconciliationCost);
		writer.println(" (Duplications: " + numberOfDuplication + ", Transfers: " + numberOfTransfer + ", Losses: "
				+ numberOfLoss + ")");
		writer.println("Total number of optimal solutions: " + getNumberOfSolutions());
		writer.println("Total number of candidates for gene birth: " + optimalMappings.get(geneTree.getRoot()).size());
		writer.close();
	}
	public void modifySolutionWithLeafCountHeuristic() {

		ArrayList<Node> allPreOrderInternalGeneNodes = geneTree.getAllPreOrderInternalNodes();
		ArrayList<Node> speciesLeaves = speciesTree.getLeafSet();
		ArrayList<Node> geneLeaves = geneTree.getLeafSet();
		HashMap<String, Integer> inferredGeneFrequencies = new HashMap<String, Integer>();
		HashMap<String, Integer> actualGeneFrequencies = getFrequencyMapFromNodes(geneLeaves);
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node speciesNode = mapping.getSpeciesNode();

			if (geneNode.getEvent().equals("Duplication")) {
			//	speciesNode.incrementDuplicationCount();
			} else if (geneNode.getEvent().equals("Transfer")) {

				speciesNode.incrementTransferCount();
				Node transferRecipient = mapping.getAdditiveTransferRecipient();
				transferRecipient.incrementTransferRecipientCount();
			}
		}
		updateMappingAndDuplicationCounts();
		initializeInferredFrequencies(inferredGeneFrequencies, speciesLeaves);
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node childCausingTransfer = null;
			if (geneNode.getOption() == 2 || geneNode.getOption() == 15)
				childCausingTransfer = geneNode.getRight();
			else
				childCausingTransfer = geneNode.getLeft();
			if (geneNode.getEvent().equals("Transfer")) {
				Node additiveRecipient = mapping.getAdditiveTransferRecipient();
				
				if (currentHeuristic(additiveRecipient, mapping, actualGeneFrequencies, inferredGeneFrequencies)) {
					geneNode.setEvent("Replacing Transfer");
				} else {
					geneNode.setEvent("Additive Transfer");

				}
			}
		}

	}

	public void modifySolutionWithMappingCountHeuristic() {

		ArrayList<Node> allPreOrderInternalGeneNodes = geneTree.getAllPreOrderInternalNodes();
		ArrayList<Node> speciesLeaves = speciesTree.getLeafSet();
		ArrayList<Node> geneLeaves = geneTree.getLeafSet();
		HashMap<String, Integer> inferredGeneFrequencies = new HashMap<String, Integer>();
		HashMap<String, Integer> actualGeneFrequencies = getFrequencyMapFromNodes(geneLeaves);
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node speciesNode = mapping.getSpeciesNode();

			if (geneNode.getEvent().equals("Duplication")) {
			//	speciesNode.incrementDuplicationCount();
			} else if (geneNode.getEvent().equals("Transfer")) {

				speciesNode.incrementTransferCount();
				Node transferRecipient = mapping.getAdditiveTransferRecipient();
				transferRecipient.incrementTransferRecipientCount();
			}
		}
		
		updateMappingAndDuplicationCounts();
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node childCausingTransfer = null;
			if (geneNode.getOption() == 2 || geneNode.getOption() == 15)
				childCausingTransfer = geneNode.getRight();
			else
				childCausingTransfer = geneNode.getLeft();
			if (geneNode.getEvent().equals("Transfer")) {
				Node additiveRecipient = mapping.getAdditiveTransferRecipient();
				/*if(additiveByGlobalMappingCount(additiveRecipient,geneNode)) {
					geneNode.setEvent("Additive Transfer");
				}
				else {
					geneNode.setEvent("Replacing Transfer");
				}*/
				
				if (additiveRecipient.additiveByMappingCount(geneNode)) {
					geneNode.setEvent("Additive Transfer");
				} else {
					geneNode.setEvent("Replacing Transfer");

				}
			}
		}


}
	private boolean additiveByGlobalMappingCount(Node additiveRecipient, GeneNodeRooted geneNodeMapped) {
		ArrayList<Node> allPreOrderInternalGeneNodes = geneTree.getAllPreOrderInternalNodes();
		int mappingCount = 0;
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node speciesNode = mapping.getSpeciesNode();
			if(speciesNode.equals(additiveRecipient)) {
				if(!geneNode.isDescendantOrEqual(geneNodeMapped) || !geneNodeMapped.isDescendantOrEqual(geneNode)) {
					mappingCount++;
				}
			}
		}
		if(mappingCount>0)
			return true;
		else
			return false;
	}

	public void modifySolutionWithLossHeuristic() {

		ArrayList<Node> allPreOrderInternalGeneNodes = geneTree.getAllPreOrderInternalNodes();
		ArrayList<Node> speciesLeaves = speciesTree.getLeafSet();
		HashMap<String, Integer> inferredGeneFrequencies = new HashMap<String, Integer>();
		
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node speciesNode = mapping.getSpeciesNode();

			if (geneNode.getEvent().equals("Duplication")) {
			//	speciesNode.incrementDuplicationCount();
			} else if (geneNode.getEvent().equals("Transfer")) {

				speciesNode.incrementTransferCount();
				Node transferRecipient = mapping.getAdditiveTransferRecipient();
				transferRecipient.incrementTransferRecipientCount();
			}
		}
		updateMappingAndDuplicationCounts();
		initializeInferredFrequencies(inferredGeneFrequencies, speciesLeaves);
		for (Node node : allPreOrderInternalGeneNodes) {
			GeneNodeRooted geneNode = (GeneNodeRooted) node;
			RootedBinaryTreeMapping mapping = geneNode.getMapping();
			Node childCausingTransfer = null;
			if (geneNode.getOption() == 2 || geneNode.getOption() == 15)
				childCausingTransfer = geneNode.getRight();
			else
				childCausingTransfer = geneNode.getLeft();
			if (geneNode.getEvent().equals("Transfer")) {
				Node additiveRecipient = mapping.getAdditiveTransferRecipient();
				
				if (additiveRecipient.anyAncestorBranchHavingLosses()) {
					geneNode.setEvent("Replacing Transfer");
				} else {
					geneNode.setEvent("Additive Transfer");

				}
			}
		}

	}



	private void updateMappingAndDuplicationCounts() {
	//	ArrayList<Node> postOrderInternalGeneNodes = geneTree.getAllPostOrderInternalNodes();
		
		for(Node geneNode: optimalMappings.keySet()) {
			ArrayList<RootedBinaryTreeMapping>candidateMappings = optimalMappings.get(geneNode);
			//sortMappingsBySpeciesNodeMapping(candidateMappings);
			for(RootedBinaryTreeMapping mapping:candidateMappings) {
				mapping.updateMappingAndDuplicationCounts();
			}
			
		}
	}

	

	private void initializeInferredFrequencies(HashMap<String, Integer> inferredGeneFrequencies,
			ArrayList<Node> speciesLeaves) {
		GeneNodeRooted geneRoot = (GeneNodeRooted) geneTree.getRoot();
		Node rootImage = geneRoot.getMapping().getSpeciesNode();

		for (Node leaf : speciesLeaves) {

			int inferredCount = leaf.getInferredFrequencyByCurrentHeuristic(rootImage);
			inferredGeneFrequencies.put(leaf.getName(), inferredCount);
		}
	}

	private boolean currentHeuristic(Node additiveRecipient, RootedBinaryTreeMapping mapping,
			HashMap<String, Integer> actualGeneFrequencies, HashMap<String, Integer> inferredGeneFrequencies) {
		ArrayList<Node> leavesUnderTransferRecipient = additiveRecipient.getAllLeavesUnderSubtree();

		for (Node leaf : leavesUnderTransferRecipient) {

			int inferredCount = inferredGeneFrequencies.get(leaf.getName());

			int actualCount = 0;
			if (actualGeneFrequencies.containsKey(leaf.getName()))
				actualCount = actualGeneFrequencies.get(leaf.getName());
			// System.out.println(leaf+ " i: "+inferredCount+ " a: "+ actualCount);
			if (inferredCount <= actualCount) {
				return false;
			}
		}
		for (Node leaf : leavesUnderTransferRecipient) {

			int presentCount = inferredGeneFrequencies.get(leaf.getName());
			inferredGeneFrequencies.put(leaf.getName(), presentCount - 1);
		}
		return true;
	}

	


	

	private Node getActualParentOfTransfer(Node additiveRecipient, ArrayList<Node> leavesCausingTransfer) {
		ArrayList<Node> leavesUnderAdditiveRecipient = additiveRecipient.getAllLeavesUnderSubtree();
		ArrayList<Node> leavesInCommon = new ArrayList<>();
		for (Node node : leavesUnderAdditiveRecipient) {
			for (Node leaf : leavesCausingTransfer) {
				if (leaf.getName().equals(node.getName())) {
					leavesInCommon.add(node);
				}
			}
		}
		Node LCA = null;
		if (leavesInCommon.size() > 0)
			LCA = leavesInCommon.get(0);
		for (int i = 1; i < leavesInCommon.size(); i++) {
			LCA = LCA.getLCA(leavesInCommon.get(i));
		}
		return LCA;
	}

	

	

	


	private HashMap<String, Integer> getFrequencyMapFromNodes(ArrayList<Node> leafSet) {
		HashMap<String, Integer> nodeFrequencyMap = new HashMap<String, Integer>();
		for (Node node : leafSet) {
			if (nodeFrequencyMap.containsKey(node.getName())) {
				int presentCount = nodeFrequencyMap.get(node.getName());
				nodeFrequencyMap.put(node.getName(), presentCount + 1);
			} else {
				nodeFrequencyMap.put(node.getName(), 1);
			}
		}
		return nodeFrequencyMap;
	}

	// DTL basic algorithm with replacement heuristic is performed here
	public void performDTLReplacementHeuristic() {

		// get all nodes of gene tree and species tree by post order
		ArrayList<Node> allPostOrderSpeciesNodes = speciesTree.getAllNodesByPostOrder();
		ArrayList<Node> allPostOrderGeneNodes = geneTree.getAllNodesByPostOrder();

		// get size of gene tree and species tree
		int geneTreeSize = allPostOrderGeneNodes.size();
		int speciesTreeSize = allPostOrderSpeciesNodes.size();

		// get leaf nodes of gene tree and species tree
		ArrayList<Node> leafGeneNodes = geneTree.getLeafSet();
		ArrayList<Node> leafSpeciesNodes = speciesTree.getLeafSet();

		// get all internal nodes of gene tree by post-order
		ArrayList<Node> internalGeneNodes = geneTree.getAllPostOrderInternalNodes();

		// get all internal nodes of species tree by pre-order
		ArrayList<Node> internalPreOrderSpeciesNodes = speciesTree.getAllPreOrderInternalNodes();

		// create a mapping from each gene node to every species node

		initializeAllMappings(geneTreeSize, speciesTreeSize, allPostOrderGeneNodes, allPostOrderSpeciesNodes);

		initializeLeafMappings(leafGeneNodes, leafSpeciesNodes);

		// for each internal gene node
		for (Node geneNode : internalGeneNodes) {

			// get left child g' and right child g" of gene node
			Node leftGeneChild = geneNode.getLeft();
			Node rightGeneChild = geneNode.getRight();

			// assume the mapping cost to be infinity first
			int minimumMappingCost = Integer.MAX_VALUE;

			// Initialize an optimal mapping for this gene node
			ArrayList<RootedBinaryTreeMapping> optimalMappingsForGeneNode = new ArrayList<>();

			// for each species node
			for (Node speciesNode : allPostOrderSpeciesNodes) {

				// let the tentative mapping be (g,s)
				RootedBinaryTreeMapping tentativeMapping = getMappingFromNodes(geneNode, speciesNode);

				// let the left child to species node mapping be (g', s)
				RootedBinaryTreeMapping leftToParentMapping = getMappingFromNodes(leftGeneChild, speciesNode);

				// let the right child to species node mapping be (g", s)
				RootedBinaryTreeMapping rightToParentMapping = getMappingFromNodes(rightGeneChild, speciesNode);

				// if s is a leaf node
				if (speciesNode.isLeaf()) {

					handleIntToLeafMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild, tentativeMapping,
							leftToParentMapping, rightToParentMapping);
					setOptionsForIntToLeafMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild,
							tentativeMapping);
				} else {
					Node leftSpeciesChild = speciesNode.getLeft();
					Node rightSpeciesChild = speciesNode.getRight();
					handleIntToIntMapping(tentativeMapping, geneNode, speciesNode, leftToParentMapping,
							rightToParentMapping);
					setOptionsForIntToIntMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild, leftSpeciesChild,
							rightSpeciesChild, tentativeMapping);
				}
				// if this mapping has smaller cost then let it be present optimal mapping
				if (tentativeMapping.getCost() < minimumMappingCost) {
					minimumMappingCost = tentativeMapping.getCost();
					optimalMappingsForGeneNode.clear();
					optimalMappingsForGeneNode.add(tentativeMapping);
				} else if (tentativeMapping.getCost() == minimumMappingCost) {
					optimalMappingsForGeneNode.add(tentativeMapping);
				}
			}
			// add the optimal mapping to optimal mapping list
			optimalMappings.put(geneNode, optimalMappingsForGeneNode);

			processOutCalculationWithHeuristic(internalPreOrderSpeciesNodes, geneNode);
		}
		// the final reconciliation cost should be the cost of mapping for root of gene
		// tree
		RootedBinaryTreeMapping rootMap = optimalMappings.get(geneTree.getRoot()).get(0);
		reconciliationCost = rootMap.getCost();
	}

	private void processOutCalculationWithHeuristic(ArrayList<Node> internalPreOrderSpeciesNodes, Node geneNode) {
		// for each species node in pre order
		for (Node speciesNode : internalPreOrderSpeciesNodes) {

			// let s' = left(s) and s" = right(s)
			Node leftSpeciesChild = speciesNode.getLeft();
			Node rightSpeciesChild = speciesNode.getRight();

			// get the mappings (g, s), (g, s'), (g, s")
			RootedBinaryTreeMapping tentativeMapping = getMappingFromNodes(geneNode, speciesNode);
			RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
			RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

			// out(g, s') = { out(g, s), inAlt(g, s")}
			int parentToLeftOutCost = Math.min(tentativeMapping.getOutCost(), parentToRightMapping.getInAltCost());
			parentToLeftMapping.setOutCost(parentToLeftOutCost);
			if (speciesNode.getLeftBranch().getNumberOfLoss() > 0 && parentToLeftOutCost < Integer.MAX_VALUE / 10
					&& parentToLeftOutCost > 0) {
				int newCost = parentToLeftOutCost - LOSS - ADDITIVE_TRANSFER + REPLACEMENT_TRANSFER;
				if (newCost < 0)
					newCost = 0;
				parentToRightMapping.setOutCost(newCost);
			}
			// out(g, s") = { out(g, s), inAlt(g, s')}
			int parentToRightOutCost = Math.min(tentativeMapping.getOutCost(), parentToLeftMapping.getInAltCost());

			parentToRightMapping.setOutCost(parentToRightOutCost);
			if (speciesNode.getRightBranch().getNumberOfLoss() > 0 && parentToRightOutCost < Integer.MAX_VALUE / 10
					&& parentToRightOutCost > 0) {
				int newCost = parentToRightOutCost - LOSS - ADDITIVE_TRANSFER + REPLACEMENT_TRANSFER;
				if (newCost < 0)
					newCost = 0;
				parentToLeftMapping.setOutCost(newCost);
			}
			updateOutSpeciesNode(tentativeMapping, parentToLeftMapping, parentToRightMapping);

		}
	}




	private String getNumberOfSolutions() {
		BigInteger totalNumberOfSolution = new BigInteger("0");
		for (RootedBinaryTreeMapping rootMapping : optimalMappings.get(geneTree.getRoot())) {
			// System.out.println(rootMapping+ ": number:
			// "+rootMapping.getNumberOfSolutions());
			// System.out.println(rootMapping.getOptions());
			// System.out.println(rootMapping.getOptionWeights());
			totalNumberOfSolution = totalNumberOfSolution.add(rootMapping.getNumberOfSolutions());
		}
		return totalNumberOfSolution.toString();
	}

	private void sampleRandomSolutionWithLossTracking() {
		ArrayList<Node> allPreOrderGeneNodes = geneTree.getAllNodesByPreOrder();
		Random random = new Random();
		for (Node geneNode : allPreOrderGeneNodes) {
			RootedBinaryTreeMapping optimalMapping = null;
			ArrayList<RootedBinaryTreeMapping> optimalMappingsOfGeneNode = optimalMappings.get(geneNode);
			GeneNodeRooted geneNodeRooted = (GeneNodeRooted) geneNode;
			if (geneNode.isRoot()) {

				optimalMapping = optimalMappingsOfGeneNode.get(random.nextInt(optimalMappingsOfGeneNode.size()));

			} else {
				GeneNodeRooted parentNode = (GeneNodeRooted) geneNodeRooted.getParent();

				RootedBinaryTreeMapping parentMapping = parentNode.getMapping();
				Node parentImage = parentMapping.getSpeciesNode();
				Node rightOfParentImage = parentMapping.getSpeciesNode().getRight();
				Node leftOfParentImage = parentMapping.getSpeciesNode().getLeft();

				int parentOption = parentNode.getOption();
				if (parentOption == 1) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 2 || parentOption == 15) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getInMappingFromPair(geneNode, parentImage);
						// parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());
					} else {
						optimalMapping = getOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 3 || parentOption == 16) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getInMappingFromPair(geneNode, parentImage);
						// parentMapping.updateInTransferLossesWithBranch(optimalMapping.getSpeciesNode());

					} else {
						optimalMapping = getOutMappingFromPair(geneNode, parentImage);
						parentNode.getMapping().setAdditiveTransferRecipient(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 4) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getInMappingFromPair(geneNode, leftOfParentImage);

					} else {

						optimalMapping = getInMappingFromPair(geneNode, rightOfParentImage);

					}
					// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 5) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getInMappingFromPair(geneNode, rightOfParentImage);
					} else {
						optimalMapping = getInMappingFromPair(geneNode, leftOfParentImage);
					}
					// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					// System.out.println( "Childnode: "+geneNode +" "+parentMapping + " losses: " +
					// parentMapping.getLosses());
				} else if (parentOption == 6) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						// parentMapping.addLoss(parentImage.getLeft());
						// parentMapping.addLostBranch(parentImage.getLeftBranch());
					} else {
						optimalMapping = getInMappingFromPair(geneNode, rightOfParentImage);// inSearchWithLossTracking(optimalMappingsOfGeneNode,
																							// rightOfParentImage);
						// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 7) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						// parentMapping.addLoss(parentImage.getRight());
						// parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getInMappingFromPair(geneNode, leftOfParentImage);
						// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
					}
				} else if (parentOption == 8) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						// parentMapping.addLoss(parentImage.getLeft());
						// parentMapping.addLostBranch(parentImage.getLeftBranch());

					} else {
						optimalMapping = getInMappingFromPair(geneNode, rightOfParentImage);
						// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 9) {
					if (parentNode.getRight().equals(geneNodeRooted)) {
						optimalMapping = getMappingFromNodes(geneNode, parentImage);
						// parentMapping.addLoss(parentImage.getRight());
						// parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getInMappingFromPair(geneNode, leftOfParentImage);
						// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

					}
				} else if (parentOption == 10) {
					optimalMapping = getMappingFromNodes(geneNode, parentImage);
				} else if (parentOption == 11) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getInMappingFromPair(geneNode, leftOfParentImage);
						// parentMapping.addLoss(parentImage.getRight());
						// parentMapping.addLostBranch(parentImage.getRightBranch());
					} else {
						optimalMapping = getInMappingFromPair(geneNode, rightOfParentImage);
						// parentMapping.addLoss(parentImage.getLeft());
						// parentMapping.addLostBranch(parentImage.getLeftBranch());
					}
					// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 12) {
					if (parentNode.getLeft().equals(geneNodeRooted)) {
						optimalMapping = getInMappingFromPair(geneNode, rightOfParentImage);
						// parentMapping.addLoss(parentImage.getRight());
						// parentMapping.addLostBranch(parentImage.getRightBranch());

					} else {
						optimalMapping = getInMappingFromPair(geneNode, leftOfParentImage);
						// parentMapping.addLoss(parentImage.getLeft());
						// parentMapping.addLostBranch(parentImage.getLeftBranch());

					}
					// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				} else if (parentOption == 13) {
					optimalMapping = getInMappingFromPair(geneNode, leftOfParentImage);
					// parentMapping.addLoss(parentImage.getRight());
					// parentMapping.addLostBranch(parentImage.getRightBranch());
					// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());
				} else if (parentOption == 14) {
					optimalMapping = getInMappingFromPair(geneNode, rightOfParentImage);
					// parentMapping.addLoss(parentImage.getLeft());
					// parentMapping.addLostBranch(parentImage.getLeftBranch());
					// parentMapping.updateLossesWithBranch(optimalMapping.getSpeciesNode());

				}
			}

			geneNodeRooted.setMappedSpeciesNode(optimalMapping.getSpeciesNode());
			geneNodeRooted
					.setOption(optimalMapping.getOptions().get(random.nextInt(optimalMapping.getOptions().size())));
			geneNodeRooted.setMapping(optimalMapping);
			geneNodeRooted.setEventByOption();

		}

	}

	private RootedBinaryTreeMapping getInMappingFromPair(Node geneNode, Node speciesNode) {
		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);
		Node mappedNode = mapping.getInSpeciesNodeRandomly();
		return getMappingFromNodes(geneNode, mappedNode);
	}

	private RootedBinaryTreeMapping getOutMappingFromPair(Node geneNode, Node speciesNode) {
		RootedBinaryTreeMapping mapping = getMappingFromNodes(geneNode, speciesNode);
		Node mappedNode = mapping.getOutSpeciesNodeRandomly();
		return getMappingFromNodes(geneNode, mappedNode);
	}

	public RootedBinaryTreeMapping getMappingFromNodes(Node geneNode, Node speciesNode) {
		// System.out.println(geneNode+ " "+speciesNode);
		return allPossibleMappings[geneNode.getPostOrderNumber()][speciesNode.getPostOrderNumber()];
	}

	// DTL basic algorithm is performed here
	public void performDTLReconciliation() {

		// get all nodes of gene tree and species tree by post order
		ArrayList<Node> allPostOrderSpeciesNodes = speciesTree.getAllNodesByPostOrder();
		ArrayList<Node> allPostOrderGeneNodes = geneTree.getAllNodesByPostOrder();

		// get size of gene tree and species tree
		int geneTreeSize = allPostOrderGeneNodes.size();
		int speciesTreeSize = allPostOrderSpeciesNodes.size();

		// get leaf nodes of gene tree and species tree
		ArrayList<Node> leafGeneNodes = geneTree.getLeafSet();
		ArrayList<Node> leafSpeciesNodes = speciesTree.getLeafSet();

		// get all internal nodes of gene tree by post-order
		ArrayList<Node> internalGeneNodes = geneTree.getAllPostOrderInternalNodes();

		// get all internal nodes of species tree by pre-order
		ArrayList<Node> internalPreOrderSpeciesNodes = speciesTree.getAllPreOrderInternalNodes();

		// create a mapping from each gene node to every species node
		allPossibleMappings = new RootedBinaryTreeMapping[geneTreeSize][speciesTreeSize];
		initializeAllMappings(geneTreeSize, speciesTreeSize, allPostOrderGeneNodes, allPostOrderSpeciesNodes);

		initializeLeafMappings(leafGeneNodes, leafSpeciesNodes);
		// get number of internal gene nodes
		int numberOfInternalGeneNodes = internalGeneNodes.size();

		// for each internal gene node
		for (int i = 0; i < numberOfInternalGeneNodes; i++) {

			// get the gene node g
			Node geneNode = internalGeneNodes.get(i);

			// get left child g' and right child g" of gene node
			Node leftGeneChild = geneNode.getLeft();
			Node rightGeneChild = geneNode.getRight();

			// assume the mapping cost to be infinity first
			int minimumMappingCost = Integer.MAX_VALUE;

			// Initialize an optimal mapping for this gene node
			ArrayList<RootedBinaryTreeMapping> optimalMappingsForSingleGeneNode = new ArrayList<>();

			// for each species node
			for (int j = 0; j < speciesTreeSize; j++) {

				// get the species node s
				Node speciesNode = allPostOrderSpeciesNodes.get(j);

				// let the tentative mapping be (g,s)
				RootedBinaryTreeMapping tentativeMapping = allPossibleMappings[geneNode.getPostOrderNumber()][j];

				// let the left child to species node mapping be (g', s)
				RootedBinaryTreeMapping leftToParentMapping = allPossibleMappings[leftGeneChild
						.getPostOrderNumber()][j];

				// let the right child to species node mapping be (g", s)
				RootedBinaryTreeMapping rightToParentMapping = allPossibleMappings[rightGeneChild
						.getPostOrderNumber()][j];

				// if s is a leaf node
				if (speciesNode.isLeaf()) {

					handleIntToLeafMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild, tentativeMapping,
							leftToParentMapping, rightToParentMapping);

				} else {
					handleIntToIntMapping(tentativeMapping, geneNode, speciesNode, leftToParentMapping,
							rightToParentMapping);

				}
				// if this mapping has smaller cost then let it be present optimal mapping
				if (tentativeMapping.getCost() < minimumMappingCost) {
					minimumMappingCost = tentativeMapping.getCost();
					optimalMappingsForSingleGeneNode.clear();
					optimalMappingsForSingleGeneNode.add(tentativeMapping);
				} else if (tentativeMapping.getCost() == minimumMappingCost) {
					optimalMappingsForSingleGeneNode.add(tentativeMapping);
				}
			}
			// add the optimal mapping to optimal mapping list
			optimalMappings.put(geneNode, optimalMappingsForSingleGeneNode);

			processOutCalculation(internalPreOrderSpeciesNodes, geneNode);

		}
		// the final reconciliation cost should be the cost of mapping for root of gene
		// tree
		reconciliationCost = optimalMappings.get(geneTree.getRoot()).get(0).getCost();
	}
	// print all the outputs (optimal mappings and reconciliation costs)

	// DTL basic algorithm is performed here
	public void performDTLReconciliationWithOptionTracking() {

		// get all nodes of gene tree and species tree by post order
		ArrayList<Node> allPostOrderSpeciesNodes = speciesTree.getAllNodesByPostOrder();
		ArrayList<Node> allPostOrderGeneNodes = geneTree.getAllNodesByPostOrder();

		// get size of gene tree and species tree
		int geneTreeSize = allPostOrderGeneNodes.size();
		int speciesTreeSize = allPostOrderSpeciesNodes.size();

		// System.out.println(speciesTreeSize);
		// get leaf nodes of gene tree and species tree
		ArrayList<Node> leafGeneNodes = geneTree.getLeafSet();
		ArrayList<Node> leafSpeciesNodes = speciesTree.getLeafSet();

		// get all internal nodes of gene tree by post-order
		ArrayList<Node> internalGeneNodes = geneTree.getAllPostOrderInternalNodes();

		// get all internal nodes of species tree by pre-order
		ArrayList<Node> internalPreOrderSpeciesNodes = speciesTree.getAllPreOrderInternalNodes();

		// create a mapping from each gene node to every species node
		initializeAllMappings(geneTreeSize, speciesTreeSize, allPostOrderGeneNodes, allPostOrderSpeciesNodes);

		initializeLeafMappings(leafGeneNodes, leafSpeciesNodes);

		// for each internal gene node
		for (Node geneNode : internalGeneNodes) {

			// get left child g' and right child g" of gene node
			Node leftGeneChild = geneNode.getLeft();
			Node rightGeneChild = geneNode.getRight();

			// assume the mapping cost to be infinity first
			int minimumMappingCost = Integer.MAX_VALUE;

			// Initialize an optimal mapping for this gene node
			ArrayList<RootedBinaryTreeMapping> optimalMappingsForGeneNode = new ArrayList<>();

			// for each species node
			for (Node speciesNode : allPostOrderSpeciesNodes) {

				// let the tentative mapping be (g,s)
				RootedBinaryTreeMapping tentativeMapping = getMappingFromNodes(geneNode, speciesNode);

				// let the left child to species node mapping be (g', s)
				RootedBinaryTreeMapping leftToParentMapping = getMappingFromNodes(leftGeneChild, speciesNode);

				// let the right child to species node mapping be (g", s)
				RootedBinaryTreeMapping rightToParentMapping = getMappingFromNodes(rightGeneChild, speciesNode);

				// if s is a leaf node
				if (speciesNode.isLeaf()) {

					handleIntToLeafMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild, tentativeMapping,
							leftToParentMapping, rightToParentMapping);
					setOptionsForIntToLeafMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild,
							tentativeMapping);
				} else {
					Node leftSpeciesChild = speciesNode.getLeft();
					Node rightSpeciesChild = speciesNode.getRight();
					handleIntToIntMapping(tentativeMapping, geneNode, speciesNode, leftToParentMapping,
							rightToParentMapping);
					setOptionsForIntToIntMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild, leftSpeciesChild,
							rightSpeciesChild, tentativeMapping);
				}

				// if this mapping has smaller cost then let it be present optimal mapping
				if (tentativeMapping.getCost() < minimumMappingCost) {
					minimumMappingCost = tentativeMapping.getCost();
					optimalMappingsForGeneNode.clear();
					optimalMappingsForGeneNode.add(tentativeMapping);
				} else if (tentativeMapping.getCost() == minimumMappingCost) {
					optimalMappingsForGeneNode.add(tentativeMapping);
				}

			}

			// add the optimal mapping to optimal mapping list

			optimalMappings.put(geneNode, optimalMappingsForGeneNode);
			// System.out.println("Before "+
			// allPossibleMappings[28][3].getOutSpeciesNode());

			processOutCalculation(internalPreOrderSpeciesNodes, geneNode);

			// }
		}
		// the final reconciliation cost should be the cost of mapping for root of gene
		// tree
		RootedBinaryTreeMapping rootMap = optimalMappings.get(geneTree.getRoot()).get(0);
		reconciliationCost = rootMap.getCost();
		// System.out.println("The reconciliation cost is: "+ reconciliationCost);
	}

	private void performAllDTLReconciliation() {

		// get all nodes of gene tree and species tree by post order
		ArrayList<Node> allPostOrderSpeciesNodes = speciesTree.getAllNodesByPostOrder();
		ArrayList<Node> allPostOrderGeneNodes = geneTree.getAllNodesByPostOrder();

		// get size of gene tree and species tree
		int geneTreeSize = allPostOrderGeneNodes.size();
		int speciesTreeSize = allPostOrderSpeciesNodes.size();

		// System.out.println(speciesTreeSize);
		// get leaf nodes of gene tree and species tree
		ArrayList<Node> leafGeneNodes = geneTree.getLeafSet();
		ArrayList<Node> leafSpeciesNodes = speciesTree.getLeafSet();

		// get all internal nodes of gene tree by post-order
		ArrayList<Node> internalGeneNodes = geneTree.getAllPostOrderInternalNodes();

		// get all internal nodes of species tree by pre-order
		ArrayList<Node> internalPreOrderSpeciesNodes = speciesTree.getAllPreOrderInternalNodes();

		// create a mapping from each gene node to every species node
		initializeAllMappings(geneTreeSize, speciesTreeSize, allPostOrderGeneNodes, allPostOrderSpeciesNodes);

		initializeLeafMappingsForMultipleOptima(leafGeneNodes, leafSpeciesNodes);

		// for each internal gene node
		for (Node geneNode : internalGeneNodes) {

			// get left child g' and right child g" of gene node
			Node leftGeneChild = geneNode.getLeft();
			Node rightGeneChild = geneNode.getRight();

			// assume the mapping cost to be infinity first
			int minimumMappingCost = Integer.MAX_VALUE;

			// Initialize an optimal mapping for this gene node
			ArrayList<RootedBinaryTreeMapping> optimalMappingsForGeneNode = new ArrayList<>();

			// for each species node
			for (Node speciesNode : allPostOrderSpeciesNodes) {

				// let the tentative mapping be (g,s)
				RootedBinaryTreeMapping tentativeMapping = getMappingFromNodes(geneNode, speciesNode);

				// let the left child to species node mapping be (g', s)
				RootedBinaryTreeMapping leftToParentMapping = getMappingFromNodes(leftGeneChild, speciesNode);

				// let the right child to species node mapping be (g", s)
				RootedBinaryTreeMapping rightToParentMapping = getMappingFromNodes(rightGeneChild, speciesNode);

				// if s is a leaf node
				if (speciesNode.isLeaf()) {

					handleIntToLeafMapping(geneNode, speciesNode, leftGeneChild, rightGeneChild, tentativeMapping,
							leftToParentMapping, rightToParentMapping);
					setOptionsForIntToLeafMappingForMultipleOptima(geneNode, speciesNode, leftGeneChild, rightGeneChild,
							tentativeMapping);
				} else {
					Node leftSpeciesChild = speciesNode.getLeft();
					Node rightSpeciesChild = speciesNode.getRight();
					handleIntToIntMappingForMultipleOptima(tentativeMapping, geneNode, speciesNode, leftToParentMapping,
							rightToParentMapping);
					setOptionsForIntToIntForMultipleOptima(geneNode, speciesNode, leftGeneChild, rightGeneChild,
							leftSpeciesChild, rightSpeciesChild, tentativeMapping);
				}

				tentativeMapping.sumUpNumberOfOptions();

				// if this mapping has smaller cost then let it be present optimal mapping
				if (tentativeMapping.getCost() < minimumMappingCost) {
					minimumMappingCost = tentativeMapping.getCost();
					optimalMappingsForGeneNode.clear();
					optimalMappingsForGeneNode.add(tentativeMapping);
				} else if (tentativeMapping.getCost() == minimumMappingCost) {
					optimalMappingsForGeneNode.add(tentativeMapping);
				}

			}
			processOutCalculationForMultipleOptima(internalPreOrderSpeciesNodes, geneNode);

			// add the optimal mapping to optimal mapping list

			optimalMappings.put(geneNode, optimalMappingsForGeneNode);

		}
		// the final reconciliation cost should be the cost of mapping for root of gene
		// tree
		RootedBinaryTreeMapping rootMap = optimalMappings.get(geneTree.getRoot()).get(0);
		reconciliationCost = rootMap.getCost();
	}

	private void handleIntToIntMappingForMultipleOptima(RootedBinaryTreeMapping tentativeMapping, Node geneNode,
			Node speciesNode, RootedBinaryTreeMapping leftToParentMapping,
			RootedBinaryTreeMapping rightToParentMapping) {

		// get left child s' and right child s" of species node
		Node leftSpeciesChild = speciesNode.getLeft();
		Node rightSpeciesChild = speciesNode.getRight();
		// get left child g' and right child g" of gene node
		Node leftGeneChild = geneNode.getLeft();
		Node rightGeneChild = geneNode.getRight();
		// let the gene node to left child of species node mapping be (g, s')
		RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

		// let the gene node to right child of species node mapping be (g, s")
		RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
		// let the left child of gene node to left child of species node mapping be (g',
		// s')
		RootedBinaryTreeMapping leftToLeftMapping = getMappingFromNodes(leftGeneChild, leftSpeciesChild);

		// let the left child of gene node to right child of species node mapping be
		// (g', s")
		RootedBinaryTreeMapping leftToRightMapping = getMappingFromNodes(leftGeneChild, rightSpeciesChild);

		// let the right child of gene node to left child of species node mapping be
		// (g", s')
		RootedBinaryTreeMapping rightToLeftMapping = getMappingFromNodes(rightGeneChild, leftSpeciesChild);

		// let the right child of gene node to right child of species node mapping be
		// (g", s")
		RootedBinaryTreeMapping rightToRightMapping = getMappingFromNodes(rightGeneChild, rightSpeciesChild);

		// speciationCost = min{in(g',s') + in (g",s"), in(g",s') + in (g', s")}
		int speciationCost = Math.min(leftToLeftMapping.getInCost() + rightToRightMapping.getInCost(),
				leftToRightMapping.getInCost() + rightToLeftMapping.getInCost());
		tentativeMapping.setSpeciationCost(speciationCost);

		// save all duplication option costs in an array

		int[] duplicationOptions = getDuplicationOptions(leftToParentMapping, rightToParentMapping, leftToLeftMapping,
				leftToRightMapping, rightToLeftMapping, rightToRightMapping);

		// Set duplicationCost = DUPLICATION + min(duplicationOptions)
		int duplicationCost = DUPLICATION + getMinimum(duplicationOptions);
		tentativeMapping.setDuplicationCost(duplicationCost);

		// s is not root
		if (!speciesNode.isRoot()) {

			// transferCost(g,s) = TRANSFER + min(in(g',s) + out(g",s), in(g",s) + out(g's))
			int transferCost = TRANSFER + Math.min(leftToParentMapping.getInCost() + rightToParentMapping.getOutCost(),
					leftToParentMapping.getOutCost() + rightToParentMapping.getInCost());
			tentativeMapping.setTransferCost(transferCost);
		}

		// c(g,s) = min(speciationCost, duplicationCost, transferCost)
		tentativeMapping.setCost(getMinimumValue(tentativeMapping.getSpeciationCost(),
				tentativeMapping.getDuplicationCost(), tentativeMapping.getTransferCost()));

		int previousInCost = tentativeMapping.getInCost();
		// in(g,s) = min{c(g,s), in(g,s') + LOSS, in(g, s") + LOSS}
		tentativeMapping.setInCost(getMinimumValue(tentativeMapping.getCost(), parentToLeftMapping.getInCost() + LOSS,
				parentToRightMapping.getInCost() + LOSS));

		updateInNodeForMultipleOptima(tentativeMapping, previousInCost);

		int previousInAltCost = tentativeMapping.getInAltCost();
		// inAlt(g,s) = min{c(g,s), inAlt(g,s'), inAlt(g, s")}
		tentativeMapping.setInAltCost(getMinimumValue(tentativeMapping.getCost(), parentToLeftMapping.getInAltCost(),
				parentToRightMapping.getInAltCost()));

		updateInAltNodeForMultipleOptima(tentativeMapping, leftSpeciesChild, rightSpeciesChild, previousInAltCost);

	}

	private void updateInAltNodeForMultipleOptima(RootedBinaryTreeMapping tentativeMapping, Node leftSpeciesChild,
			Node rightSpeciesChild, int previousCost) {

		if (tentativeMapping.getInAltCost() < previousCost) {
			tentativeMapping.getInAltSpeciesNodes().clear();
		}
		Node geneNode = tentativeMapping.getGeneNode();
		// get left child s' and right child s" of species node

		// let the gene node to left child of species node mapping be (g, s')
		RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

		// let the gene node to right child of species node mapping be (g, s")
		RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
		if (tentativeMapping.getInAltCost() == parentToRightMapping.getInAltCost()) {
			for (Node inAltNode : parentToRightMapping.getInAltSpeciesNodes()) {
				if (!tentativeMapping.getInAltSpeciesNodes().contains(inAltNode))
					tentativeMapping.addInAltSpeciesNode(inAltNode);
			}
		}
		if (tentativeMapping.getInAltCost() == parentToLeftMapping.getInAltCost()) {
			for (Node inAltNode : parentToLeftMapping.getInAltSpeciesNodes())

			{
				if (!tentativeMapping.getInAltSpeciesNodes().contains(inAltNode))
					tentativeMapping.addInAltSpeciesNode(inAltNode);
			}
		}
		if (tentativeMapping.getInAltCost() == tentativeMapping.getCost()) {
			if (!tentativeMapping.getInAltSpeciesNodes().contains(tentativeMapping.getSpeciesNode()))
				tentativeMapping.addInAltSpeciesNode(tentativeMapping.getSpeciesNode());
		}

	}

	private void updateInNodeForMultipleOptima(RootedBinaryTreeMapping tentativeMapping, int previousInCost) {

		if (tentativeMapping.getInCost() < previousInCost) {
			tentativeMapping.getInSpeciesNodes().clear();
		}
		Node speciesNode = tentativeMapping.getSpeciesNode();
		Node geneNode = tentativeMapping.getGeneNode();
		// get left child s' and right child s" of species node
		Node leftSpeciesChild = speciesNode.getLeft();
		Node rightSpeciesChild = speciesNode.getRight();
		// let the gene node to left child of species node mapping be (g, s')
		RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

		// let the gene node to right child of species node mapping be (g, s")
		RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
		if (tentativeMapping.getInCost() == parentToLeftMapping.getInCost() + LOSS) {
			for (Node inNode : parentToLeftMapping.getInSpeciesNodes()) {
				if (!tentativeMapping.getInSpeciesNodes().contains(inNode)) {
					tentativeMapping.addInSpeciesNode(inNode);
				}
			}

		}
		if (tentativeMapping.getInCost() == parentToRightMapping.getInCost() + LOSS) {
			for (Node inNode : parentToRightMapping.getInSpeciesNodes()) {
				if (!tentativeMapping.getInSpeciesNodes().contains(inNode)) {
					tentativeMapping.addInSpeciesNode(inNode);
				}
			}
		}
		if (tentativeMapping.getInCost() == tentativeMapping.getCost()) {
			if (!tentativeMapping.getInSpeciesNodes().contains(speciesNode))
				tentativeMapping.addInSpeciesNode(speciesNode);
		}

	}

	private void processOutCalculationForMultipleOptima(ArrayList<Node> internalPreOrderSpeciesNodes, Node geneNode) {

		// for each species node in pre order
		for (Node speciesNode : internalPreOrderSpeciesNodes) {

			// let s' = left(s) and s" = right(s)
			Node leftSpeciesChild = speciesNode.getLeft();
			Node rightSpeciesChild = speciesNode.getRight();

			// get the mappings (g, s), (g, s'), (g, s")
			RootedBinaryTreeMapping tentativeMapping = getMappingFromNodes(geneNode, speciesNode);
			RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
			RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

			// int prevParentToLeftOutCost = parentToLeftMapping.getOutCost();
			// out(g, s') = { out(g, s), inAlt(g, s")}
			parentToLeftMapping
					.setOutCost(Math.min(tentativeMapping.getOutCost(), parentToRightMapping.getInAltCost()));
			// if (parentToLeftMapping.getOutCost() < prevParentToLeftOutCost) {
			parentToLeftMapping.getOutSpeciesNodes().clear();
			// }

			// int prevParentToRightOutCost = parentToRightMapping.getOutCost();
			// out(g, s") = { out(g, s), inAlt(g, s')}
			parentToRightMapping
					.setOutCost(Math.min(tentativeMapping.getOutCost(), parentToLeftMapping.getInAltCost()));
			// if (parentToRightMapping.getOutCost() < prevParentToRightOutCost) {
			parentToRightMapping.getOutSpeciesNodes().clear();
			// }
			updateOutSpeciesNodeForMultipleOptima(tentativeMapping, parentToLeftMapping, parentToRightMapping);

		}

	}

	private void updateOutSpeciesNodeForMultipleOptima(RootedBinaryTreeMapping tentativeMapping,
			RootedBinaryTreeMapping parentToLeftMapping, RootedBinaryTreeMapping parentToRightMapping) {

		if (parentToLeftMapping.getOutCost() == tentativeMapping.getOutCost()) {
			for (Node outNode : tentativeMapping.getOutSpeciesNodes()) {
				if (!parentToLeftMapping.getOutSpeciesNodes().contains(outNode))
					parentToLeftMapping.addOutSpeciesNode(outNode);
			}
		}
		if (parentToLeftMapping.getOutCost() == parentToRightMapping.getInAltCost()) {
			for (Node inAltNode : parentToRightMapping.getInAltSpeciesNodes()) {
				if (!parentToLeftMapping.getOutSpeciesNodes().contains(inAltNode))
					parentToLeftMapping.addOutSpeciesNode(inAltNode);
			}
		}
		if (parentToRightMapping.getOutCost() == tentativeMapping.getOutCost()) {
			for (Node outNode : tentativeMapping.getOutSpeciesNodes()) {
				if (!parentToRightMapping.getOutSpeciesNodes().contains(outNode))
					parentToRightMapping.addOutSpeciesNode(outNode);
			}
		}
		if (parentToRightMapping.getOutCost() == parentToLeftMapping.getInAltCost()) {
			for (Node inAltNode : parentToLeftMapping.getInAltSpeciesNodes()) {
				if (!parentToRightMapping.getOutSpeciesNodes().contains(inAltNode))
					parentToRightMapping.addOutSpeciesNode(inAltNode);
			}
		}

	}

	private void setOptionsForIntToIntForMultipleOptima(Node geneNode, Node speciesNode, Node leftGeneChild,
			Node rightGeneChild, Node leftSpeciesChild, Node rightSpeciesChild,
			RootedBinaryTreeMapping tentativeMapping) {

		// let the left child to species node mapping be (g', s)
		RootedBinaryTreeMapping leftToParentMapping = getMappingFromNodes(leftGeneChild, speciesNode);

		// let the right child to species node mapping be (g", s)
		RootedBinaryTreeMapping rightToParentMapping = getMappingFromNodes(rightGeneChild, speciesNode);

		// let the left child of gene node to left child of species node mapping be (g',
		// s')
		RootedBinaryTreeMapping leftToLeftMapping = getMappingFromNodes(leftGeneChild, leftSpeciesChild);

		// let the left child of gene node to right child of species node mapping be
		// (g', s")
		RootedBinaryTreeMapping leftToRightMapping = getMappingFromNodes(leftGeneChild, rightSpeciesChild);

		// let the right child of gene node to left child of species node mapping be
		// (g", s')
		RootedBinaryTreeMapping rightToLeftMapping = getMappingFromNodes(rightGeneChild, leftSpeciesChild);

		// let the right child of gene node to right child of species node mapping be
		// (g", s")
		RootedBinaryTreeMapping rightToRightMapping = getMappingFromNodes(rightGeneChild, rightSpeciesChild);

		// speciationCost = min{in(g',s') + in (g",s"), in(g",s') + in (g', s")}
		int speciationCost = Math.min(leftToLeftMapping.getInCost() + rightToRightMapping.getInCost(),
				leftToRightMapping.getInCost() + rightToLeftMapping.getInCost());
		if (tentativeMapping.getCost() >= Integer.MAX_VALUE / 10)
			return;
		if (tentativeMapping.getCost() == speciationCost) {
			if (tentativeMapping.getCost() == leftToLeftMapping.getInCost() + rightToRightMapping.getInCost()) {
				tentativeMapping.getOptions().add(4);
				tentativeMapping.addWeightByInAndIn(leftToLeftMapping, rightToRightMapping, this);
			}
			if (tentativeMapping.getCost() == leftToRightMapping.getInCost() + rightToLeftMapping.getInCost()) {
				tentativeMapping.getOptions().add(5);
				tentativeMapping.addWeightByInAndIn(leftToRightMapping, rightToLeftMapping, this);
			}
		}

		int[] duplicationOptions = getDuplicationOptions(leftToParentMapping, rightToParentMapping, leftToLeftMapping,
				leftToRightMapping, rightToLeftMapping, rightToRightMapping);

		// Set duplicationCost = DUPLICATION + min(duplicationOptions)
		int duplicationCost = DUPLICATION + getMinimum(duplicationOptions);
		if (tentativeMapping.getCost() == duplicationCost) {
			if (tentativeMapping.getCost() - DUPLICATION == leftToParentMapping.getCost()
					+ rightToRightMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(6);
				tentativeMapping.addWeightByNormalAndIn(leftToParentMapping, rightToRightMapping, this);
			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToParentMapping.getCost()
					+ rightToLeftMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(7);
				tentativeMapping.addWeightByNormalAndIn(leftToParentMapping, rightToLeftMapping, this);

			}
			if (tentativeMapping.getCost() - DUPLICATION == rightToParentMapping.getCost()
					+ leftToRightMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(8);
				tentativeMapping.addWeightByNormalAndIn(rightToParentMapping, leftToRightMapping, this);
			}
			if (tentativeMapping.getCost() - DUPLICATION == rightToParentMapping.getCost()
					+ leftToLeftMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(9);
				tentativeMapping.addWeightByNormalAndIn(rightToParentMapping, leftToLeftMapping, this);

			}
			if (tentativeMapping.getCost() - DUPLICATION == rightToParentMapping.getCost()
					+ leftToParentMapping.getCost()) {
				tentativeMapping.getOptions().add(10);
				tentativeMapping.addWeightByNormalAndNormal(rightToParentMapping, leftToParentMapping);

			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToLeftMapping.getInCost()
					+ rightToRightMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(11);
				tentativeMapping.addWeightByInAndIn(leftToLeftMapping, rightToRightMapping, this);
			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToRightMapping.getInCost()
					+ rightToLeftMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(12);
				tentativeMapping.addWeightByInAndIn(leftToRightMapping, rightToLeftMapping, this);

			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToLeftMapping.getInCost()
					+ rightToLeftMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(13);
				tentativeMapping.addWeightByInAndIn(leftToLeftMapping, rightToLeftMapping, this);

			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToRightMapping.getInCost()
					+ rightToRightMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(14);
				tentativeMapping.addWeightByInAndIn(leftToRightMapping, rightToRightMapping, this);

			}
		}
		if (!speciesNode.isRoot()) {

			// transferCost(g,s) = TRANSFER + min(in(g',s) + out(g",s), in(g",s) + out(g's))
			int transferCost = TRANSFER + Math.min(leftToParentMapping.getInCost() + rightToParentMapping.getOutCost(),
					leftToParentMapping.getOutCost() + rightToParentMapping.getInCost());
			if (tentativeMapping.getCost() == transferCost) {
				if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getInCost()
						+ rightToParentMapping.getOutCost()) {
					tentativeMapping.getOptions().add(15);
					tentativeMapping.addWeightByInOut(leftToParentMapping, rightToParentMapping, this);
				}
				if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getOutCost()
						+ rightToParentMapping.getInCost()) {
					tentativeMapping.getOptions().add(16);
					tentativeMapping.addWeightByInOut(rightToParentMapping, leftToParentMapping, this);
				}
			}

		}

	}

	private void setOptionsForIntToLeafMappingForMultipleOptima(Node geneNode, Node speciesNode, Node leftGeneChild,
			Node rightGeneChild, RootedBinaryTreeMapping tentativeMapping) {

		RootedBinaryTreeMapping leftToParentMapping = getMappingFromNodes(leftGeneChild, speciesNode);

		// let the right child to species node mapping be (g", s)
		RootedBinaryTreeMapping rightToParentMapping = getMappingFromNodes(rightGeneChild, speciesNode);
		int duplicationCost = DUPLICATION + leftToParentMapping.getCost() + rightToParentMapping.getCost();

		if (tentativeMapping.getCost() == duplicationCost) {
			tentativeMapping.getOptions().add(1);
			tentativeMapping.addWeightByChildren(leftToParentMapping, rightToParentMapping);
		}
		// if s is not root
		if (!speciesNode.isRoot()) {

			// transferCost(g,s) = TRANSFER + min(in(g',s) + out(g",s), in(g",s) + out(g's))
			if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getInCost()
					+ rightToParentMapping.getOutCost()) {
				tentativeMapping.getOptions().add(2);
				tentativeMapping.addWeightByInOut(leftToParentMapping, rightToParentMapping, this);
			}
			if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getOutCost()
					+ rightToParentMapping.getInCost()) {
				tentativeMapping.getOptions().add(3);
				tentativeMapping.addWeightByInOut(rightToParentMapping, leftToParentMapping, this);
			}
		}

	}

	private void initializeLeafMappingsForMultipleOptima(ArrayList<Node> leafGeneNodes,
			ArrayList<Node> leafSpeciesNodes) {

		// for each leaf node of gene tree
		for (Node geneLeafNode : leafGeneNodes) {
			// for each leaf node of species tree
			for (Node speciesLeafNode : leafSpeciesNodes) {

				// gene node name and species node names are equal
				if (geneLeafNode.getName().equals(speciesLeafNode.getName())) {

					RootedBinaryTreeMapping mapping = getMappingFromNodes(geneLeafNode, speciesLeafNode);
					// set options to 1
					mapping.setNumberOfSolutions(new BigInteger("1"));
					// set cost and in-cost to zero
					mapping.setCost(0);
					mapping.setInCost(0);
					mapping.setInAltCost(0);

			//		speciesLeafNode.incrementMappingCount();

					ArrayList<RootedBinaryTreeMapping> optimalMappingsForSingleGeneNode = new ArrayList<>();

					optimalMappingsForSingleGeneNode.add(mapping);
					mapping.getOptions().add(0);
					mapping.getOptionWeights().add(new BigInteger("1"));
					mapping.addInSpeciesNode(speciesLeafNode);
					mapping.addInAltSpeciesNode(speciesLeafNode);
					// allPossibleMappings[leafGeneNodeIndex][ancestorSpeciesNodeIndex].getOptions().add(0);
					// add leaf to leaf match to optimal mappings
					optimalMappings.put(geneLeafNode, optimalMappingsForSingleGeneNode);
					// get all ancestors of this species node
					ArrayList<Node> ancestorSpeciesNodes = speciesLeafNode.getAllAncestors();

					// for each ancestor s' of this species node
					for (Node ancestorSpeciesNode : ancestorSpeciesNodes) {

						// get distance from species node to ancestor node
						int distance = speciesLeafNode.getDistance(ancestorSpeciesNode);

						// set inAlt(g, s') = 0 and in(g, s') = distance*LOSS
						RootedBinaryTreeMapping ancestorMapping = getMappingFromNodes(geneLeafNode,
								ancestorSpeciesNode);
						ancestorMapping.setInCost(LOSS * distance);
						ancestorMapping.setInAltCost(0);
						ancestorMapping.addInSpeciesNode(speciesLeafNode);
						ancestorMapping.addInAltSpeciesNode(speciesLeafNode);

					}

				}
			}
		}

		for (Node geneLeafNode : leafGeneNodes) {
			//System.out.println(geneLeafNode);
			//System.out.println(geneLeafNode);
			Node actualMappedNode = optimalMappings.get(geneLeafNode).get(0).getSpeciesNode();

			for (Node species : speciesTree.getAllNodesByPostOrder()) {
				if (!species.equals(actualMappedNode) && !actualMappedNode.isDescendant(species)) {
					RootedBinaryTreeMapping mapping = getMappingFromNodes(geneLeafNode, species);
					mapping.addOutSpeciesNode(actualMappedNode);
					mapping.setOutCost(0);
				}
			}
		}

	}

	private void updateInAltNodeForIntToIntMapping(RootedBinaryTreeMapping tentativeMapping, Node leftSpeciesChild,
			Node rightSpeciesChild, int previousCost) {
		if (tentativeMapping.getInAltCost() < previousCost) {
			tentativeMapping.getInAltSpeciesNodes().clear();
		}
		Node geneNode = tentativeMapping.getGeneNode();
		// get left child s' and right child s" of species node

		// let the gene node to left child of species node mapping be (g, s')
		RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

		// let the gene node to right child of species node mapping be (g, s")
		RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
		if (tentativeMapping.getInAltCost() == parentToRightMapping.getInAltCost()) {
			tentativeMapping.addInAltSpeciesNode(parentToRightMapping.getInAltSpeciesNodeRandomly());
		}
		if (tentativeMapping.getInAltCost() == parentToLeftMapping.getInAltCost()) {
			tentativeMapping.addInAltSpeciesNode(parentToLeftMapping.getInAltSpeciesNodeRandomly());
		}
		if (tentativeMapping.getInAltCost() == tentativeMapping.getCost()) {
			tentativeMapping.addInAltSpeciesNode(tentativeMapping.getSpeciesNode());
		}
	}

	private void updateInNodeForIntToIntMapping(RootedBinaryTreeMapping tentativeMapping, int previousInCost) {
		// if(tentativeMapping.getInCost()<previousInCost) {
		tentativeMapping.getInSpeciesNodes().clear();
		// }
		Node speciesNode = tentativeMapping.getSpeciesNode();
		Node geneNode = tentativeMapping.getGeneNode();
		// get left child s' and right child s" of species node
		Node leftSpeciesChild = speciesNode.getLeft();
		Node rightSpeciesChild = speciesNode.getRight();
		// let the gene node to left child of species node mapping be (g, s')
		RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

		// let the gene node to right child of species node mapping be (g, s")
		RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
		if (tentativeMapping.getInCost() == parentToLeftMapping.getInCost() + LOSS) {
			tentativeMapping.addInSpeciesNode(parentToLeftMapping.getInSpeciesNodeRandomly());
		}
		if (tentativeMapping.getInCost() == parentToRightMapping.getInCost() + LOSS) {
			tentativeMapping.addInSpeciesNode(parentToRightMapping.getInSpeciesNodeRandomly());
		}
		if (tentativeMapping.getInCost() == tentativeMapping.getCost()) {
			tentativeMapping.addInSpeciesNode(speciesNode);
		}

	}

	private void handleIntToIntMapping(RootedBinaryTreeMapping tentativeMapping, Node geneNode, Node speciesNode,
			RootedBinaryTreeMapping leftToParentMapping, RootedBinaryTreeMapping rightToParentMapping) {
		// get left child s' and right child s" of species node
		Node leftSpeciesChild = speciesNode.getLeft();
		Node rightSpeciesChild = speciesNode.getRight();
		// get left child g' and right child g" of gene node
		Node leftGeneChild = geneNode.getLeft();
		Node rightGeneChild = geneNode.getRight();
		// let the gene node to left child of species node mapping be (g, s')
		RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

		// let the gene node to right child of species node mapping be (g, s")
		RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
		// let the left child of gene node to left child of species node mapping be (g',
		// s')
		RootedBinaryTreeMapping leftToLeftMapping = getMappingFromNodes(leftGeneChild, leftSpeciesChild);

		// let the left child of gene node to right child of species node mapping be
		// (g', s")
		RootedBinaryTreeMapping leftToRightMapping = getMappingFromNodes(leftGeneChild, rightSpeciesChild);

		// let the right child of gene node to left child of species node mapping be
		// (g", s')
		RootedBinaryTreeMapping rightToLeftMapping = getMappingFromNodes(rightGeneChild, leftSpeciesChild);

		// let the right child of gene node to right child of species node mapping be
		// (g", s")
		RootedBinaryTreeMapping rightToRightMapping = getMappingFromNodes(rightGeneChild, rightSpeciesChild);

		// speciationCost = min{in(g',s') + in (g",s"), in(g",s') + in (g', s")}
		int speciationCost = Math.min(leftToLeftMapping.getInCost() + rightToRightMapping.getInCost(),
				leftToRightMapping.getInCost() + rightToLeftMapping.getInCost());
		tentativeMapping.setSpeciationCost(speciationCost);

		// save all duplication option costs in an array

		int[] duplicationOptions = getDuplicationOptions(leftToParentMapping, rightToParentMapping, leftToLeftMapping,
				leftToRightMapping, rightToLeftMapping, rightToRightMapping);

		// Set duplicationCost = DUPLICATION + min(duplicationOptions)
		int duplicationCost = DUPLICATION + getMinimum(duplicationOptions);
		tentativeMapping.setDuplicationCost(duplicationCost);

		// s is not root
		if (!speciesNode.isRoot()) {

			// transferCost(g,s) = TRANSFER + min(in(g',s) + out(g",s), in(g",s) + out(g's))
			int transferCost = TRANSFER + Math.min(leftToParentMapping.getInCost() + rightToParentMapping.getOutCost(),
					leftToParentMapping.getOutCost() + rightToParentMapping.getInCost());
			tentativeMapping.setTransferCost(transferCost);
		}

		// c(g,s) = min(speciationCost, duplicationCost, transferCost)
		tentativeMapping.setCost(getMinimumValue(tentativeMapping.getSpeciationCost(),
				tentativeMapping.getDuplicationCost(), tentativeMapping.getTransferCost()));

		int previousInCost = tentativeMapping.getInCost();
		// in(g,s) = min{c(g,s), in(g,s') + LOSS, in(g, s") + LOSS}
		tentativeMapping.setInCost(getMinimumValue(tentativeMapping.getCost(), parentToLeftMapping.getInCost() + LOSS,
				parentToRightMapping.getInCost() + LOSS));

		updateInNodeForIntToIntMapping(tentativeMapping, previousInCost);

		int previousInAltCost = tentativeMapping.getInAltCost();
		// inAlt(g,s) = min{c(g,s), inAlt(g,s'), inAlt(g, s")}
		tentativeMapping.setInAltCost(getMinimumValue(tentativeMapping.getCost(), parentToLeftMapping.getInAltCost(),
				parentToRightMapping.getInAltCost()));

		updateInAltNodeForIntToIntMapping(tentativeMapping, leftSpeciesChild, rightSpeciesChild, previousInAltCost);

	}

	private void processOutCalculation(ArrayList<Node> internalPreOrderSpeciesNodes, Node geneNode) {
		// for each species node in pre order
		for (Node speciesNode : internalPreOrderSpeciesNodes) {

			// let s' = left(s) and s" = right(s)
			Node leftSpeciesChild = speciesNode.getLeft();
			Node rightSpeciesChild = speciesNode.getRight();

			// get the mappings (g, s), (g, s'), (g, s")
			RootedBinaryTreeMapping tentativeMapping = getMappingFromNodes(geneNode, speciesNode);
			RootedBinaryTreeMapping parentToLeftMapping = getMappingFromNodes(geneNode, leftSpeciesChild);
			RootedBinaryTreeMapping parentToRightMapping = getMappingFromNodes(geneNode, rightSpeciesChild);

			int prevParentToLeftOutCost = parentToLeftMapping.getOutCost();
			// out(g, s') = { out(g, s), inAlt(g, s")}
			parentToLeftMapping
					.setOutCost(Math.min(tentativeMapping.getOutCost(), parentToRightMapping.getInAltCost()));
			if (parentToLeftMapping.getOutCost() < prevParentToLeftOutCost) {
				parentToLeftMapping.getOutSpeciesNodes().clear();
			}

			int prevParentToRightOutCost = parentToRightMapping.getOutCost();
			// out(g, s") = { out(g, s), inAlt(g, s')}
			parentToRightMapping
					.setOutCost(Math.min(tentativeMapping.getOutCost(), parentToLeftMapping.getInAltCost()));
			if (parentToRightMapping.getOutCost() < prevParentToRightOutCost) {
				parentToRightMapping.getOutSpeciesNodes().clear();
			}
			updateOutSpeciesNode(tentativeMapping, parentToLeftMapping, parentToRightMapping);

		}
	}

	private void updateOutSpeciesNode(RootedBinaryTreeMapping tentativeMapping,
			RootedBinaryTreeMapping parentToLeftMapping, RootedBinaryTreeMapping parentToRightMapping) {

		if (parentToLeftMapping.getOutCost() == tentativeMapping.getOutCost()) {
			parentToLeftMapping.addOutSpeciesNode(tentativeMapping.getOutSpeciesNodeRandomly());
		} else {
			parentToLeftMapping.addOutSpeciesNode(parentToRightMapping.getInAltSpeciesNodeRandomly());
		}
		if (parentToRightMapping.getOutCost() == tentativeMapping.getOutCost()) {
			parentToRightMapping.addOutSpeciesNode(tentativeMapping.getOutSpeciesNodeRandomly());
		} else {
			parentToRightMapping.addOutSpeciesNode(parentToLeftMapping.getInAltSpeciesNodeRandomly());
		}

	}

	private int[] getDuplicationOptions(RootedBinaryTreeMapping leftToParentMapping,
			RootedBinaryTreeMapping rightToParentMapping, RootedBinaryTreeMapping leftToLeftMapping,
			RootedBinaryTreeMapping leftToRightMapping, RootedBinaryTreeMapping rightToLeftMapping,
			RootedBinaryTreeMapping rightToRightMapping) {

		int[] duplicationOptions = { leftToParentMapping.getCost() + rightToRightMapping.getInCost() + LOSS,
				leftToParentMapping.getCost() + rightToLeftMapping.getInCost() + LOSS,
				rightToParentMapping.getCost() + leftToRightMapping.getInCost() + LOSS,
				rightToParentMapping.getCost() + leftToLeftMapping.getInCost() + LOSS,
				rightToParentMapping.getCost() + leftToParentMapping.getCost(),
				leftToLeftMapping.getInCost() + rightToRightMapping.getInCost() + 2 * LOSS,
				leftToRightMapping.getInCost() + rightToLeftMapping.getInCost() + 2 * LOSS,
				leftToLeftMapping.getInCost() + rightToLeftMapping.getInCost() + 2 * LOSS,
				leftToRightMapping.getInCost() + rightToRightMapping.getInCost() + 2 * LOSS };

		return duplicationOptions;
	}

	private void handleIntToLeafMapping(Node geneNode, Node speciesNode, Node leftGeneChild, Node rightGeneChild,
			RootedBinaryTreeMapping tentativeMapping, RootedBinaryTreeMapping leftToParentMapping,
			RootedBinaryTreeMapping rightToParentMapping) {
		// internal node to leaf node speciation cost should be infinity
		tentativeMapping.setSpeciationCost(Integer.MAX_VALUE / 10);

		// duplicationCost(g,s) = DUPLICATION + c(g',s) + c(g",s)
		int duplicationCost = DUPLICATION + leftToParentMapping.getCost() + rightToParentMapping.getCost();
		tentativeMapping.setDuplicationCost(duplicationCost);

		// if s is not root
		if (!speciesNode.isRoot()) {

			// transferCost(g,s) = TRANSFER + min(in(g',s) + out(g",s), in(g",s) + out(g's))
			int transferCost = TRANSFER + Math.min(leftToParentMapping.getInCost() + rightToParentMapping.getOutCost(),
					leftToParentMapping.getOutCost() + rightToParentMapping.getInCost());
			tentativeMapping.setTransferCost(transferCost);

		}

		// c(g,s) = min(speciationCost, duplicationCost, transferCost)
		tentativeMapping.setCost(getMinimumValue(tentativeMapping.getSpeciationCost(),
				tentativeMapping.getDuplicationCost(), tentativeMapping.getTransferCost()));

		// in(g,s) = c(g,s)
		tentativeMapping.setInCost(tentativeMapping.getCost());

		tentativeMapping.getInSpeciesNodes().clear();
		tentativeMapping.addInSpeciesNode(speciesNode);

		// inAlt(g,s) = c(g,s)
		tentativeMapping.setInAltCost(tentativeMapping.getCost());

		tentativeMapping.getInAltSpeciesNodes().clear();
		tentativeMapping.addInAltSpeciesNode(speciesNode);

	}

	private void initializeAllMappings(int geneTreeSize, int speciesTreeSize, ArrayList<Node> allPostOrderGeneNodes,
			ArrayList<Node> allPostOrderSpeciesNodes) {
		allPossibleMappings = new RootedBinaryTreeMapping[geneTreeSize][speciesTreeSize];
		for (int i = 0; i < geneTreeSize; i++) {
			Node geneNode = allPostOrderGeneNodes.get(i);
			for (int j = 0; j < speciesTreeSize; j++) {

				Node speciesNode = allPostOrderSpeciesNodes.get(j);
				allPossibleMappings[i][j] = new RootedBinaryTreeMapping(geneNode, speciesNode);

			}

		}

	}

	private void initializeLeafMappings(ArrayList<Node> leafGeneNodes, ArrayList<Node> leafSpeciesNodes) {

		// for each leaf node of gene tree
		for (Node geneLeafNode : leafGeneNodes) {
			// for each leaf node of species tree
			for (Node speciesLeafNode : leafSpeciesNodes) {

				// gene node name and species node names are equal
				if (geneLeafNode.getName().equals(speciesLeafNode.getName())) {
					RootedBinaryTreeMapping mapping = getMappingFromNodes(geneLeafNode, speciesLeafNode);
					// set cost and in-cost to zero
					mapping.setCost(0);
					mapping.setInCost(0);
					mapping.setInAltCost(0);

					ArrayList<RootedBinaryTreeMapping> optimalMappingsForSingleGeneNode = new ArrayList<>();

					optimalMappingsForSingleGeneNode.add(mapping);
					mapping.getOptions().add(0);
					mapping.addInSpeciesNode(speciesLeafNode);
					mapping.addInAltSpeciesNode(speciesLeafNode);
					// allPossibleMappings[leafGeneNodeIndex][ancestorSpeciesNodeIndex].getOptions().add(0);
					// add leaf to leaf match to optimal mappings
					optimalMappings.put(geneLeafNode, optimalMappingsForSingleGeneNode);
					// get all ancestors of this species node
					ArrayList<Node> ancestorSpeciesNodes = speciesLeafNode.getAllAncestors();

					// for each ancestor s' of this species node
					for (Node ancestorSpeciesNode : ancestorSpeciesNodes) {

						// get distance from species node to ancestor node
						int distance = speciesLeafNode.getDistance(ancestorSpeciesNode);

						// set inAlt(g, s') = 0 and in(g, s') = distance*LOSS
						RootedBinaryTreeMapping ancestorMapping = getMappingFromNodes(geneLeafNode,
								ancestorSpeciesNode);
						ancestorMapping.setInCost(LOSS * distance);
						ancestorMapping.setInAltCost(0);
						ancestorMapping.addInSpeciesNode(speciesLeafNode);
						ancestorMapping.addInAltSpeciesNode(speciesLeafNode);

					}

				}
			}
		}

		for (Node geneLeafNode : leafGeneNodes) {
			Node actualMappedNode = optimalMappings.get(geneLeafNode).get(0).getSpeciesNode();

			for (Node species : speciesTree.getAllNodesByPostOrder()) {
				if (!species.equals(actualMappedNode) && !actualMappedNode.isDescendant(species)) {
					RootedBinaryTreeMapping mapping = getMappingFromNodes(geneLeafNode, species);
					mapping.addOutSpeciesNode(actualMappedNode);
					mapping.setOutCost(0);
				}
			}
		}
	}

	private void setOptionsForIntToIntMapping(Node geneNode, Node speciesNode, Node leftGeneChild, Node rightGeneChild,
			Node leftSpeciesChild, Node rightSpeciesChild, RootedBinaryTreeMapping tentativeMapping) {

		// let the left child to species node mapping be (g', s)
		RootedBinaryTreeMapping leftToParentMapping = getMappingFromNodes(leftGeneChild, speciesNode);

		// let the right child to species node mapping be (g", s)
		RootedBinaryTreeMapping rightToParentMapping = getMappingFromNodes(rightGeneChild, speciesNode);

		// let the left child of gene node to left child of species node mapping be (g',
		// s')
		RootedBinaryTreeMapping leftToLeftMapping = getMappingFromNodes(leftGeneChild, leftSpeciesChild);

		// let the left child of gene node to right child of species node mapping be
		// (g', s")
		RootedBinaryTreeMapping leftToRightMapping = getMappingFromNodes(leftGeneChild, rightSpeciesChild);

		// let the right child of gene node to left child of species node mapping be
		// (g", s')
		RootedBinaryTreeMapping rightToLeftMapping = getMappingFromNodes(rightGeneChild, leftSpeciesChild);

		// let the right child of gene node to right child of species node mapping be
		// (g", s")
		RootedBinaryTreeMapping rightToRightMapping = getMappingFromNodes(rightGeneChild, rightSpeciesChild);

		// speciationCost = min{in(g',s') + in (g",s"), in(g",s') + in (g', s")}
		int speciationCost = Math.min(leftToLeftMapping.getInCost() + rightToRightMapping.getInCost(),
				leftToRightMapping.getInCost() + rightToLeftMapping.getInCost());
		if (tentativeMapping.getCost() >= Integer.MAX_VALUE / 10)
			return;
		if (tentativeMapping.getCost() == speciationCost) {
			if (tentativeMapping.getCost() == leftToLeftMapping.getInCost() + rightToRightMapping.getInCost()) {
				tentativeMapping.getOptions().add(4);
			}
			if (tentativeMapping.getCost() == leftToRightMapping.getInCost() + rightToLeftMapping.getInCost()) {
				tentativeMapping.getOptions().add(5);
			}
		}

		int[] duplicationOptions = getDuplicationOptions(leftToParentMapping, rightToParentMapping, leftToLeftMapping,
				leftToRightMapping, rightToLeftMapping, rightToRightMapping);

		// Set duplicationCost = DUPLICATION + min(duplicationOptions)
		int duplicationCost = DUPLICATION + getMinimum(duplicationOptions);
		if (tentativeMapping.getCost() == duplicationCost) {
			if (tentativeMapping.getCost() - DUPLICATION == leftToParentMapping.getCost()
					+ rightToRightMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(6);
			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToParentMapping.getCost()
					+ rightToLeftMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(7);
			}
			if (tentativeMapping.getCost() - DUPLICATION == rightToParentMapping.getCost()
					+ leftToRightMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(8);
			}
			if (tentativeMapping.getCost() - DUPLICATION == rightToParentMapping.getCost()
					+ leftToLeftMapping.getInCost() + LOSS) {
				tentativeMapping.getOptions().add(9);
			}
			if (tentativeMapping.getCost() - DUPLICATION == rightToParentMapping.getCost()
					+ leftToParentMapping.getCost()) {
				tentativeMapping.getOptions().add(10);
			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToLeftMapping.getInCost()
					+ rightToRightMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(11);
			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToRightMapping.getInCost()
					+ rightToLeftMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(12);
			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToLeftMapping.getInCost()
					+ rightToLeftMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(13);
			}
			if (tentativeMapping.getCost() - DUPLICATION == leftToRightMapping.getInCost()
					+ rightToRightMapping.getInCost() + 2 * LOSS) {
				tentativeMapping.getOptions().add(14);
			}
		}
		if (!speciesNode.isRoot()) {

			// transferCost(g,s) = TRANSFER + min(in(g',s) + out(g",s), in(g",s) + out(g's))
			int transferCost = TRANSFER + Math.min(leftToParentMapping.getInCost() + rightToParentMapping.getOutCost(),
					leftToParentMapping.getOutCost() + rightToParentMapping.getInCost());
			if (tentativeMapping.getCost() == transferCost) {
				if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getInCost()
						+ rightToParentMapping.getOutCost())
					tentativeMapping.getOptions().add(15);
				if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getOutCost()
						+ rightToParentMapping.getInCost())
					tentativeMapping.getOptions().add(16);
			}

		}
	}

	private void setOptionsForIntToLeafMapping(Node geneNode, Node speciesNode, Node leftGeneChild, Node rightGeneChild,
			RootedBinaryTreeMapping tentativeMapping) {
		RootedBinaryTreeMapping leftToParentMapping = getMappingFromNodes(leftGeneChild, speciesNode);

		// let the right child to species node mapping be (g", s)
		RootedBinaryTreeMapping rightToParentMapping = getMappingFromNodes(rightGeneChild, speciesNode);
		int duplicationCost = DUPLICATION + leftToParentMapping.getCost() + rightToParentMapping.getCost();

		if (tentativeMapping.getCost() == duplicationCost)
			tentativeMapping.getOptions().add(1);
		// if s is not root
		if (!speciesNode.isRoot()) {

			// transferCost(g,s) = TRANSFER + min(in(g',s) + out(g",s), in(g",s) + out(g's))

			if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getInCost()
					+ rightToParentMapping.getOutCost())
				tentativeMapping.getOptions().add(2);
			if (tentativeMapping.getCost() == TRANSFER + leftToParentMapping.getOutCost()
					+ rightToParentMapping.getInCost())
				tentativeMapping.getOptions().add(3);
		}

	}

	// get minimum value of an array
	private int getMinimum(int[] array) {
		int minimum = Integer.MAX_VALUE;
		for (int value : array) {
			if (value < minimum)
				minimum = value;
		}
		return minimum;
	}

	// get minimum value of three integers
	private int getMinimumValue(int a, int b, int c) {

		return Math.min(a, Math.min(b, c));
	}

	// perform some pre-processing of gene tree and species tree by pre-order,
	// post-order and in-order
	public void performPreprocessings() {
		// speciesTree.computeAllPreOrderNumbers();
		speciesTree.computeAllPostOrderNumbers();
		speciesTree.computeAllInOrderNumbers();
		speciesTree.computeLeafSet();
		// speciesTree.computeAllPreOrderInternalNodes();
		speciesTree.computeAllPostOrderInternalNodes();
		// geneTree.computeAllPreOrderNumbers();
		geneTree.computeAllPostOrderNumbers();
		geneTree.computeAllInOrderNumbers();
		geneTree.computeLeafSet();
		geneTree.computeAllPreOrderInternalNodes();
		geneTree.computeAllPostOrderInternalNodes();

		speciesTree.setAllBranches();
		geneTree.setAllBranches();
	}
}
