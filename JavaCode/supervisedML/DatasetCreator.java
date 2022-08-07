package supervisedML;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Scanner;

import dtl_algo.DTLReconciliation;
import node.Node;
import simulatedData.CreateSimulatedInput;
import utility.FileUtility;

public class DatasetCreator {

	public static void main(String[] args) throws IOException {
		
		//createRealDatasetInput();
		//combineOptRootWithSpeciesTree();
		createSupervisedMLInput();
	
	}

	private static void combineOptRootWithSpeciesTree() throws FileNotFoundException {
		Scanner speciesTreeScanner = new Scanner(new File("speciesLina.txt"));
		String speciesTree = speciesTreeScanner.nextLine();
		
		File folder = new File("OPtrootedTrees");
		File[] listOfFiles = folder.listFiles();
		int i = 0;
		for (i =0; i<listOfFiles.length;i++) {
			File optrootedFile = listOfFiles[i];
			Scanner geneTreeScanner = new Scanner(optrootedFile);
			String geneTree = geneTreeScanner.nextLine();
			String inputData = speciesTree + '\n' + geneTree;
			PrintWriter writer = new PrintWriter(new File("RealData2/input"+i+".txt"));
			writer.print(inputData);
			geneTreeScanner.close();
			writer.close();
			i++;
			
		}
		speciesTreeScanner.close();
	}

	private static void combineOptRootOutputs() throws FileNotFoundException {
		File folder = new File("RangerInputRealData2");
		File[] listOfFiles = folder.listFiles();
		for (File rangerInputFile : listOfFiles) {
			Scanner rangerScanner = new Scanner(rangerInputFile);
			Scanner optrootScanner = new Scanner(new File("OptRootOutputs/" + rangerInputFile.getName()));
			String speciesTreeString = rangerScanner.nextLine();
			for (int i = 0; i < 3; i++)
				optrootScanner.nextLine();
			String geneTreeString = optrootScanner.nextLine();
			String combinedOutput = speciesTreeString+"\n"+geneTreeString;
			PrintWriter writer = new PrintWriter(new File("OptRootedRealData2/"+rangerInputFile.getName()));
			writer.print(combinedOutput);
			writer.close();
			optrootScanner.close();
			rangerScanner.close();
		}
	}

	private static void runOptRoot() {
		try {
			// Command to create an external process
			String command = "OptRoot.win -i RealData2/input0.txt -o optRootOutput/output.txt";

			// Running the above command
			Runtime run = Runtime.getRuntime();
			Process proc = run.exec(command);
		}

		catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static void createSpecialCrossValidationSet() throws IOException {

		for (int i = 0; i < 5; i++) {
			String fileName = "all_training";
			fileName = fileName + i + ".txt";
			// createHeader(fileName);
			String datasetList[] = { "species200_" };
			File folder = new File(datasetList[0] + i);
			File[] listOfFiles = folder.listFiles();
			for (File file : listOfFiles) {
				String folderNames[] = { "00", "01", "02", "10", "11", "12", "20", "21", "22" };
				for (String folderName : folderNames) {

					if (file.getName().startsWith("dtl_trees_" + folderName)) {
						fileName = getTestFileFromFolder(folderName);
						if (i == 0) {
							createHeader(fileName);
						}
						File[] filesInSubFolder = file.listFiles();
						for (File inputFile : filesInSubFolder) {
							if (inputFile.getName().startsWith("dtl_input")) {
								String trueLabelFileName = inputFile.getPath().replaceAll("dtl", "gene")
										.replace("input", "tree") + ".pruned.guest2host";
								String arguments[] = { "-i", inputFile.getPath(), "-o", fileName, "-real",
										trueLabelFileName };
								DTLReconciliation.main(arguments);
							} else {
								inputFile.delete();
							}
						}

					} else if (file.getName().startsWith("gene_trees_")) {
						File[] filesInSubFolder = file.listFiles();
						for (File inputFile : filesInSubFolder) {
							if (!inputFile.getName().contains(".pruned.guest2host")) {
								inputFile.delete();
							}
						}
					}
				}
			}

		}
	}

	private static void createCrossValidationDataset() throws IOException {

		String datasetList[] = { "first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "nineth",
				"tenth" };
		for (int i = 0; i < 10; i++) {
			String fileName = "all_training";
			fileName = fileName + i + ".txt";
			createHeader(fileName);
			for (String datasetNumber : datasetList) {

				if (datasetNumber.equals(datasetList[i])) {
					createTestDataset(i);
					continue;
				}
				File folder = new File(datasetNumber + "_dataset");
				File[] listOfFiles = folder.listFiles();
				for (File file : listOfFiles) {
					String folderNames[] = { "00", "01", "02", "10", "11", "12", "20", "21", "22" };
					for (String folderName : folderNames) {

						if (file.getName().startsWith("dtl_trees_" + folderName)) {
							// fileName = getTestFileFromFolder(folderName);
							// createHeader(fileName);
							File[] filesInSubFolder = file.listFiles();
							for (File inputFile : filesInSubFolder) {
								if (inputFile.getName().startsWith("dtl_input")) {
									String trueLabelFileName = inputFile.getPath().replaceAll("dtl", "gene")
											.replace("input", "tree") + ".pruned.guest2host";
									String arguments[] = { "-i", inputFile.getPath(), "-o", fileName, "-real",
											trueLabelFileName };
									DTLReconciliation.main(arguments);
								} else {
									inputFile.delete();
								}
							}

						} else if (file.getName().startsWith("gene_trees_")) {
							File[] filesInSubFolder = file.listFiles();
							for (File inputFile : filesInSubFolder) {
								if (!inputFile.getName().contains(".pruned.guest2host")) {
									inputFile.delete();
								}
							}
						}
					}
				}

			}
		}
	}

	private static void createTestDataset(int i) throws IOException {
		String fileName = "all_training.txt";
		new File("test" + i).mkdirs();
		String datasetList[] = { "first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "nineth",
				"tenth" };
		String datasetNumber = datasetList[i];
		File folder = new File(datasetNumber + "_dataset");
		File[] listOfFiles = folder.listFiles();
		for (File file : listOfFiles) {
			String folderNames[] = { "00", "01", "02", "10", "11", "12", "20", "21", "22" };
			for (String folderName : folderNames) {

				if (file.getName().startsWith("dtl_trees_" + folderName)) {
					fileName = getTestFileFromFolder(folderName);
					fileName = "test" + i + "/" + fileName;
					createHeader(fileName);
					File[] filesInSubFolder = file.listFiles();
					for (File inputFile : filesInSubFolder) {
						if (inputFile.getName().startsWith("dtl_input")) {
							String trueLabelFileName = inputFile.getPath().replaceAll("dtl", "gene").replace("input",
									"tree") + ".pruned.guest2host";
							String arguments[] = { "-i", inputFile.getPath(), "-o", fileName, "-real",
									trueLabelFileName };
							DTLReconciliation.main(arguments);
						} else {
							inputFile.delete();
						}
					}

				} else if (file.getName().startsWith("gene_trees_")) {
					File[] filesInSubFolder = file.listFiles();
					for (File inputFile : filesInSubFolder) {
						if (!inputFile.getName().contains(".pruned.guest2host")) {
							inputFile.delete();
						}
					}
				}
			}
		}

	}

	private static void getAverageDTLRate(String fileName) throws FileNotFoundException {
		double totalLeaf = 0;
		double totalDuplication = 0;
		double totalTransfer = 0;
		Scanner scanner = new Scanner(new File(fileName));
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			String[] lineParts = line.split(" ");
			int numberOfLeaves = Integer.parseInt(lineParts[0]);
			int numberOfDuplications = Integer.parseInt(lineParts[1]);
			int numberOfTransfers = Integer.parseInt(lineParts[2]);

			totalLeaf = totalLeaf + numberOfLeaves;
			totalDuplication += numberOfDuplications;
			totalTransfer += numberOfTransfers;

		}
		double averageLeaf = totalLeaf / 900.0;
		double averageDuplications = totalDuplication / 900.0;
		double averageTransfers = totalTransfer / 900.0;
		System.out.println("Number of leaves on average: " + averageLeaf);
		System.out.println("Number of duplications on average: " + averageDuplications);
		System.out.println("Number of transfers on average: " + averageTransfers);
		System.out.println("Duplication ratio: " + averageDuplications / averageLeaf * 100.0 + " %");
		System.out.println("Transfer ratio: " + averageTransfers / averageLeaf * 100.0 + " %");

		scanner.close();
	}

	private static void getAverageGeneTreeSize(String file) throws FileNotFoundException {
		double total = 0;
		Scanner scanner = new Scanner(new File(file));
		while (scanner.hasNextLine()) {
			int value = Integer.parseInt(scanner.nextLine());
			total = total + value;
		}
		double average = total / 900.0;
		System.out.println(average);
		scanner.close();
	}

	private static void createRealDatasetInput() throws IOException {
		File folder = new File("realData2");
		File[] listOfFiles = folder.listFiles();
		File optrootFolder = new File("OPtrootedTrees");
		File[] optRootFiles = optrootFolder.listFiles();
		createRealHeader("realData2Across.txt");
		for (int i = 0; i<listOfFiles.length;i++) {
			File inputFile = listOfFiles[i];
			String optrootFileName = optRootFiles[i].getName();
			//System.out.println(inputFile.getName());
			String arguments[] = { "-i", inputFile.getPath(), "-o", "realData2Across.txt", "-real", optrootFileName };
			DTLReconciliation.main(arguments);
/*
			DTLReconciliation reconciliation = new DTLReconciliation();
			FileUtility.processInput(inputFile.getPath(), reconciliation);
			Node speciesRoot = reconciliation.getSpeciesTree().getRoot();
			Node geneRoot = reconciliation.getGeneTree().getRoot();
			String speciesTreeString = CreateSimulatedInput.treeToNewick(speciesRoot);
			String geneTreeString = CreateSimulatedInput.treeToNewick(geneRoot);
			String rangerInputFormat = speciesTreeString + "\n" + geneTreeString;
			PrintWriter writer = new PrintWriter(new File("RangerInputRealData2/" + inputFile.getName()));
			writer.print(rangerInputFormat);
			writer.close();
*/
			// speciesRoot
		}
	}

	public static void createRealDataInputFiles() throws IOException {
		String speciesTreeText = new String(Files.readAllBytes(Paths.get("RAXML_Trees/species.newick")));
		File folder = new File("RAXML_Trees/3GeneTrees");
		File[] listOfFiles = folder.listFiles();
		int i = 0;
		for (File inputGeneTreeFile : listOfFiles) {

			String geneTreeText = new String(Files.readAllBytes(Paths.get(inputGeneTreeFile.getPath())));
			if (geneTreeText.isEmpty() || geneTreeText.length() < 2)
				continue;
			String inputText = speciesTreeText + geneTreeText;
			PrintWriter writer = new PrintWriter(new File("RealData2/input" + i + ".txt"));
			i++;
			writer.print(inputText);
			writer.close();
			// if(i%1000==0)System.out.println(i);

		}
		System.out.println(speciesTreeText);
	}

	private static void createSupervisedMLInput() throws IOException {
		String fileName = "species_25_0.txt";
		//createHeader(fileName);
		String datasetList[] = {"first","second","third",
				 "fourth", "fifth", "tenth", "seventh", "eighth","nineth","tenth"
				};
		for (String datasetNumber : datasetList) {
			File folder = new File(datasetNumber+"_dataset" );
			File[] listOfFiles = folder.listFiles();
			for (File file : listOfFiles) {
				String folderNames[] = { "00", "01", "02", "10", "11", "12", "20", "21", "22" };
				for (String folderName : folderNames) {

					if (file.getName().startsWith("dtl_trees_" + folderName)) {
						fileName = getTestFileFromFolder(folderName);
						createHeader(fileName);
						File[] filesInSubFolder = file.listFiles();
						for (File inputFile : filesInSubFolder) {
							if (inputFile.getName().startsWith("dtl_input")) {
								String trueLabelFileName = inputFile.getPath().replaceAll("dtl", "gene")
										.replace("input", "tree") + ".pruned.guest2host";
								String arguments[] = { "-i", inputFile.getPath(), "-o", fileName, "-real",
										trueLabelFileName };
								DTLReconciliation.main(arguments);
							} else {
								inputFile.delete();
							}
						}

					} else if (file.getName().startsWith("gene_trees_")) {
						File[] filesInSubFolder = file.listFiles();
						for (File inputFile : filesInSubFolder) {
							if (!inputFile.getName().contains(".pruned.guest2host")) {
								inputFile.delete();
							}
						}
					}
				}
			}

		}
	}

	private static void outputCurrentHeuristicAccuracies() throws FileNotFoundException {

		String folderNames[] = { "00", "01", "02", "10", "11", "12", "20", "21", "22" };
		for (String folderName : folderNames) {
			String fileName = getTestFileFromFolder(folderName);
			System.out.println(fileName);
			Scanner scanner = new Scanner(new File(fileName));
			int truePositive = 0;
			int trueNegative = 0;
			int falsePositive = 0;
			int falseNegative = 0;
			scanner.nextLine();
			while (scanner.hasNextLine()) {
				String transferRowData = scanner.nextLine();
				String[] parts = transferRowData.split(",");
				int trueLabel = Integer.parseInt(parts[parts.length - 1]);
				int currentHeuristicLabel = Integer.parseInt(parts[parts.length - 2]);
				if (trueLabel == 1 && currentHeuristicLabel == 1) {
					truePositive++;
				}
				if (trueLabel == 1 && currentHeuristicLabel != 1) {
					falseNegative++;
				}
				if (trueLabel == 0 && currentHeuristicLabel == 0) {
					trueNegative++;
				}
				if (trueLabel == 0 && currentHeuristicLabel != 0) {
					falsePositive++;
				}
			}
			double additiveAccuracy = truePositive / (truePositive + falseNegative + .00001);
			double replacingAccuracy = trueNegative / (falsePositive + trueNegative + 0.00001);
			System.out.println("Current heuristic Additive accuracy: " + additiveAccuracy);
			System.out.println("Current heuristic Replacing accuracy: " + replacingAccuracy + "\n");
			scanner.close();
		}
	}

	private static void outputLossHeuristicAccuracies() throws FileNotFoundException {
		String folderNames[] = { "00", "01", "02", "10", "11", "12", "20", "21", "22" };
		for (String folderName : folderNames) {
			String fileName = getTestFileFromFolder(folderName);
			System.out.println(fileName);
			Scanner scanner = new Scanner(new File(fileName));
			int truePositive = 0;
			int trueNegative = 0;
			int falsePositive = 0;
			int falseNegative = 0;
			scanner.nextLine();

			while (scanner.hasNextLine()) {
				String transferRowData = scanner.nextLine();
				String[] parts = transferRowData.split(",");
				int trueLabel = Integer.parseInt(parts[parts.length - 1]);
				int lossHeuristicLabel = Integer.parseInt(parts[parts.length - 3]);
				if (trueLabel == 1 && lossHeuristicLabel == 1) {
					truePositive++;
				}
				if (trueLabel == 1 && lossHeuristicLabel != 1) {
					falseNegative++;
				}
				if (trueLabel == 0 && lossHeuristicLabel == 0) {
					trueNegative++;
				}
				if (trueLabel == 0 && lossHeuristicLabel != 0) {
					falsePositive++;
				}
			}
			double additiveAccuracy = truePositive / (truePositive + falseNegative + .00001);
			double replacingAccuracy = trueNegative / (falsePositive + trueNegative + 0.00001);
			System.out.println("Loss heuristic Additive accuracy: " + additiveAccuracy);
			System.out.println("Loss heuristic Replacing accuracy: " + replacingAccuracy + "\n");
			scanner.close();
		}
	}

	private static void createRealHeader(String fileName) throws IOException {
		FileWriter writer = new FileWriter(fileName, true);
		//original version
		
		writer.write("id,height,mappingCountHeuristic,lostGene2Heuristic,lostGene1Heuristic,lostGene3Heuristic,"
				+ "geneFrequencyHeuristic,distanceOfTransfer,"  + "trueLabel\r\n");
		 
		//writer.write("geneFamilyName,donor,recipient,mappingCountHeuristic,lossHeuristic,currentHeuristicLabel\n");
		writer.close();
	}
	private static void createHeader(String fileName) throws IOException {
		FileWriter writer = new FileWriter(fileName, true);
		//original version		
		writer.write("id,height,mappingCountHeuristic,lostGene2Heuristic,lostGene1Heuristic,lostGene3Heuristic,"
				+ "geneFrequencyHeuristic,"  + "trueLabel\r\n");
		
		 
		//writer.write("geneFamilyName,donor,recipient,mappingCountHeuristic,lossHeuristic,currentHeuristicLabel\n");
		writer.close();
	}

	private static String getTestFileFromFolder(String folderName) {
		if (folderName.equals("00")) {
			return "low_additive_test.txt";
		}
		if (folderName.equals("01")) {
			return "low_mixed_test.txt";
		}
		if (folderName.equals("02")) {
			return "low_replacing_test.txt";
		}
		if (folderName.equals("10")) {
			return "medium_additive_test.txt";
		}
		if (folderName.equals("11")) {
			return "medium_mixed_test.txt";
		}
		if (folderName.equals("12")) {
			return "medium_replacing_test.txt";
		}
		if (folderName.equals("20")) {
			return "high_additive_test.txt";
		}
		if (folderName.equals("21")) {
			return "high_mixed_test.txt";
		}
		if (folderName.equals("22")) {
			return "high_replacing_test.txt";
		}
		return null;
	}

	private static void createDistanceMatrixWithTrueLabel() throws IOException {
		File folder = new File("new_dataset");
		File[] listOfFiles = folder.listFiles();
		for (File file : listOfFiles) {
			if (file.getName().startsWith("dtl_trees")) {
				new File("mlInput/" + file.getName()).mkdirs();
				File[] filesInSubFolder = file.listFiles();
				for (File inputFile : filesInSubFolder) {
					if (inputFile.getName().startsWith("dtl_input")) {
						String trueLabelFileName = inputFile.getPath().replaceAll("dtl", "gene").replace("input",
								"tree") + ".pruned.guest2host";
						String arguments[] = { "-i", inputFile.getPath(), "-o",
								"mlInput/" + file.getName() + "/" + inputFile.getName(), "-real", trueLabelFileName };
						DTLReconciliation.main(arguments);
					}
				}

			}
		}
	}

	private static void runGreedyCurrentHeuristic() throws IOException {
		File folder = new File("new_dataset");
		File[] listOfFiles = folder.listFiles();
		for (File file : listOfFiles) {
			if (file.getName().startsWith("dtl_trees")) {
				new File("mlInput/" + file.getName()).mkdirs();
				File[] filesInSubFolder = file.listFiles();
				for (File inputFile : filesInSubFolder) {
					if (inputFile.getName().startsWith("dtl_input")) {
						String arguments[] = { "-i", inputFile.getPath(), "-o",
								"mlInput/" + file.getName() + "/" + inputFile.getName() };
						DTLReconciliation.main(arguments);
					}
				}

			}
		}
	}

}
