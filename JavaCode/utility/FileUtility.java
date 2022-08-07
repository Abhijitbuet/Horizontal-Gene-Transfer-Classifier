package utility;

import java.io.*;
import java.nio.charset.Charset;
import java.nio.file.*;
import java.util.*;

import dtl_algo.DTLReconciliation;
import node.*;

import tree.*;

public class FileUtility {
	// create species tree from input string (species tree part)
	public static SpeciesTreeRooted getSpeciesTreeFromFile(String treeString) throws IOException {

		treeString = treeString.replaceAll("\\s", "");

		// System.out.println(treeString);
		int length = treeString.length();
		Stack<SpeciesNodeRooted> stack = new Stack<SpeciesNodeRooted>();
		for (int i = 0; i < length; i++) {
			char presentCharacter = treeString.charAt(i);
			if (presentCharacter == '(' || presentCharacter == ',' || presentCharacter == ';')
				continue;
			else if (presentCharacter == ')') {
				SpeciesNodeRooted speciesRight = stack.pop();
				SpeciesNodeRooted speciesLeft = stack.pop();

				int nextStopCharIndex = getNextStopCharIndex(treeString, i);
				SpeciesNodeRooted parent = new SpeciesNodeRooted(speciesLeft.getName() + "$" + speciesRight.getName());
				
				if (nextStopCharIndex != i + 1 && nextStopCharIndex != -1) {
					String newParentName = (treeString.substring(i + 1, nextStopCharIndex));
					parent.setName(refineNameFromNewick(newParentName));
					i = nextStopCharIndex - 1;
				}
				parent.setLeft(speciesLeft);
				parent.setRight(speciesRight);
				// System.out.println(parent.getName());
				stack.push(parent);
			} else {
				String species = "";
				while (treeString.charAt(i) != ',' && treeString.charAt(i) != ')') {
					species += treeString.charAt(i);
					i++;
				}
				if (treeString.charAt(i) == ')')
					i--;
				String fullName = species;
				
					
				species = refineNameFromNewick(species);
				SpeciesNodeRooted speciesNode = new SpeciesNodeRooted(species);
				if(fullName.contains("_")) {
					speciesNode.setRealName(fullName);
				}
				
				stack.push(speciesNode);
				//System.out.println(speciesNode.getName());
			}
		}
		SpeciesNodeRooted root = stack.pop();
		SpeciesTreeRooted speciesTree = new SpeciesTreeRooted(root);
		return speciesTree;
	}

	private static String refineNameFromNewick(String species) {
		if(species.contains(":"))
			species = species.split(":")[0];
		
		species = species.replaceAll("_", "");
		species = species.replaceAll("-", "");
		return species;
	}

	private static int getNextStopCharIndex(String treeString, int i) {
		int nextCommaIndex = treeString.indexOf(',', i + 1);
		int nextSemiColonIndex = treeString.length();
		int nextBracketIndex = treeString.indexOf(')', i + 1);

		if (nextBracketIndex == -1)
			nextBracketIndex = Integer.MAX_VALUE;
		if (nextCommaIndex == -1)
			nextCommaIndex = Integer.MAX_VALUE;
		if (nextSemiColonIndex == -1)
			nextSemiColonIndex = Integer.MAX_VALUE;
		// System.out.println(nextSemiColonIndex);
		return Math.min(nextCommaIndex, Math.min(nextBracketIndex, nextSemiColonIndex));
	}

	public static String alignRangerSpeciesTreeFormat(Node node) {

		String output = "";
		String leftPart = "";
		String rightPart = "";
		if (node.getLeft().isLeaf()) {
			leftPart = node.getLeft().getName();
		} else {
			leftPart = alignRangerSpeciesTreeFormat(node.getLeft());
		}
		if (node.getRight().isLeaf()) {
			rightPart = node.getRight().getName();
		} else {
			if (node.getLeft().isLeaf())
				rightPart = alignRangerSpeciesTreeFormat(node.getRight());
			else
				rightPart = alignRangerSpeciesTreeFormat(node.getRight());

		}
		output = "(" + leftPart + "," + rightPart + ")" + node.getName();
		if (node.isRoot()) {

			output += ";";
		}

		return output;

	}

	public static void processInput(String fileName, DTLReconciliation reconciliation) throws IOException {
		//byte[] encoded = Files.readAllBytes(Paths.get(fileName));
		//String s = new String(encoded, Charset.defaultCharset());
		
		Scanner scanner = new Scanner(new File(fileName));
		String s = "";
		while(scanner.hasNextLine()) {
			s = s + scanner.nextLine();
		}
		scanner.close();
		s = s.replaceAll("\\s", "");
		String[] splitPartsOfInput = s.split(";");
		if (splitPartsOfInput.length < 2) {
			reconciliation.setInputError(true);
			System.out.println("input error");
			System.out.println(fileName);
			return;
		}
		String speciesTreeString = splitPartsOfInput[0];
		String geneTreeString = splitPartsOfInput[1];

		if (geneTreeString.length() == 0 || speciesTreeString.length() == 0) {
			reconciliation.setInputError(true);
			System.out.println("input error");
			System.out.println(fileName);
			return;
		}
		SpeciesTreeRooted speciesTree = getSpeciesTreeFromFile(speciesTreeString);
		GeneTreeRooted geneTree = getGeneTreeFromFile(geneTreeString);

		setPreOrderNames(speciesTree);
		String rangerSpeciesTreeFormat = alignRangerSpeciesTreeFormat(speciesTree.getRoot());
		speciesTree.setNewickOutputFormat(rangerSpeciesTreeFormat);
		// System.out.println(rangerSpeciesTreeFormat);

		geneTree.computeAllPreOrderNumbers();
		String rangerGeneTreeFormat = alignRangerGeneTreeFormat(geneTree.getRoot());
		geneTree.setNewickOutputFormat(rangerGeneTreeFormat);
		// System.out.println(rangerGeneTreeFormat);

		reconciliation.setSpeciesTree(speciesTree);
		reconciliation.setGeneTree(geneTree);
	}

	public static void setPreOrderNames(RootedTree speciesTree) {
		speciesTree.computeAllPreOrderNumbers();
		ArrayList<Node> allInternalPreOrderNodes = speciesTree.computeAllPreOrderInternalNodes();

		int index = 1;
		for (Node node : allInternalPreOrderNodes) {
			if (node.getName().contains("$"))
				node.setName("n" + index);
			index++;
		}
		// speciesTree.getAllNodesByPreOrder().clear();;
	}

	public static String alignRangerGeneTreeFormat(Node node) {

		String name = "m" + (node.getPreOrderNumber() + 1);
		node.setName(name);
		String output = "";
		String leftPart = "";
		String rightPart = "";
		if (node.getLeft().isLeaf()) {
			leftPart = node.getLeft().getName();
		} else {
			leftPart = alignRangerGeneTreeFormat(node.getLeft());
		}
		if (node.getRight().isLeaf()) {
			rightPart = node.getRight().getName();
		} else {
			if (node.getLeft().isLeaf())
				rightPart = alignRangerGeneTreeFormat(node.getRight());
			else
				rightPart = alignRangerGeneTreeFormat(node.getRight());

		}
		output = "(" + leftPart + "," + rightPart + ")" + name;
		if (node.isRoot()) {

			output += ";";
		}

		return output;
	}

//create gene tree from input string (gene tree part)
	private static GeneTreeRooted getGeneTreeFromFile(String geneTreeString) {
		geneTreeString = geneTreeString.replaceAll("\\s", "");
		geneTreeString = geneTreeString + ";";
		boolean flag = false;
		int length = geneTreeString.length();
		Stack<GeneNodeRooted> stack = new Stack<GeneNodeRooted>();
		String gene = "";
		GeneNodeRooted parent = null;
		for (int i = 0; i < length; i++) {
			
			char presentCharacter = geneTreeString.charAt(i);
			if (presentCharacter == '(' || presentCharacter == ',') {
				flag = false;
				continue;
			} else if (presentCharacter == ')') {
				GeneNodeRooted geneRight = stack.pop();
				GeneNodeRooted geneLeft = stack.pop();
				 parent = new GeneNodeRooted(geneLeft.getName() + "$" + geneRight.getName());

				//parent.setInternalNodeName(gene);
				parent.setLeft(geneLeft);
				parent.setRight(geneRight);

				stack.push(parent);
				flag = true;
			} else {
				gene = "";
				while (geneTreeString.charAt(i) != ',' && geneTreeString.charAt(i) != ')'
						&& geneTreeString.charAt(i) != ';') {
					gene += geneTreeString.charAt(i);
					i++;
				}

				if (geneTreeString.charAt(i) == ')')
					i--;
				if (gene.contains(".")) {
					String[] partsOfGeneName = gene.split(".");
					String geneName = partsOfGeneName[0];
					String realName = partsOfGeneName[1];
					GeneNodeRooted geneNode = new GeneNodeRooted(geneName);
					geneNode.setRealName(realName);
					if (!flag)
						stack.push(geneNode);
					else
						{
						parent.setInternalNodeName(gene);
						flag = false;
						}
				} else if (gene.contains("_")) {
					String[] partsOfGeneName = gene.split("_");
					String geneName = partsOfGeneName[0];
					String realName = partsOfGeneName[1];
					GeneNodeRooted geneNode = new GeneNodeRooted(geneName);
					geneNode.setRealName(realName);
					if (!flag)
						stack.push(geneNode);
					else
						{
						flag = false;
						parent.setInternalNodeName(gene);
						
						}
				} else {
					GeneNodeRooted geneNode = new GeneNodeRooted(gene);
					if (!flag)
						stack.push(geneNode);
					else
						{
						parent.setInternalNodeName(gene);						
						flag = false;
						}
				}
			}
		}
		if (stack.isEmpty())
			System.out.println("Tree: " + geneTreeString.length());
		GeneNodeRooted root = stack.pop();
		GeneTreeRooted geneTree = new GeneTreeRooted(root);
		return geneTree;
	}
}
