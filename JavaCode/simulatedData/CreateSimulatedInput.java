package simulatedData;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Scanner;
import java.util.Stack;

import node.GeneNodeRooted;
import node.Node;
import node.SpeciesNodeRooted;
import tree.GeneTreeRooted;
import tree.RootedTree;
import tree.SpeciesTreeRooted;
import utility.FileUtility;

public class CreateSimulatedInput {

	public static void main(String[] args) throws IOException {
		writeOriginalCostToFiles(1,2,3,2);
		for(int i = 1; i<=100; i++) {
			String path = "simulated/species"+i+"/";
			/*Scanner speciesTreeScanner = new Scanner(new File(path+"species.pruned.tree"));
			
			String speciesTreeString = speciesTreeScanner.nextLine();
			RootedTree speciesTree = getSpeciesTreeFromSimulatedData(speciesTreeString);
			speciesTreeScanner.close();*/
			for(int j = 1; j<=10; j++) {
			/*	Scanner geneTreeScanner = new Scanner(new File(path+ "gene"+j+".pruned.tree"));
				String geneTreeString = geneTreeScanner.nextLine();
				
				RootedTree geneTree = getGeneTreeFromSimulatedData(geneTreeString);
				PrintWriter writer = new PrintWriter(new File(path+"species"+i+"gene"+j+".newick"));
				writer.println(treeToNewick(speciesTree.getRoot()));
				writer.println(treeToNewick(geneTree.getRoot()));	
				writer.close();
				geneTreeScanner.close();*/
			
			}
		}
	}
	private static void writeOriginalCostToFiles(int loss, int duplicationCost, int additiveTransfer, int replacementTransfer ) throws FileNotFoundException {
		String outputFolderPath = "simulated\\costInfo "+ loss+ " "+duplicationCost+ " "+ additiveTransfer+ " "+ replacementTransfer+"\\";
		for(int i = 1; i<=100; i++) {
			String path = "simulated/species"+i+"/";
			
			
			for(int j = 1; j<=10; j++) {
				Scanner infoScanner = new Scanner(new File(path+"gene"+j+".pruned.info"));
				for(int k=0;k<6;k++) {
					infoScanner.nextLine();
				}
				String duplicationLine = infoScanner.nextLine();
				int numOfDuplication = Integer.parseInt(duplicationLine.split("No. of duplications:	")[1]);
				
				String lossesLine = infoScanner.nextLine();
				int numOfLosses = Integer.parseInt(lossesLine.split("No. of losses:	")[1]);
				
				String replacingLossesLine = infoScanner.nextLine();
				numOfLosses += Integer.parseInt(replacingLossesLine.split("No. of replacing losses:	")[1]);
				
				String additiveTransfersLine = infoScanner.nextLine();
				int numOfAdditiveTransfers= Integer.parseInt(additiveTransfersLine.split("No. of additive transfers:	")[1]);
				
				String replacingTransfersLine = infoScanner.nextLine();
				int numOfReplacingTransfers= Integer.parseInt(replacingTransfersLine.split("No. of replacing transfers:	")[1]);
				
				int reconciliationCost = numOfLosses*loss + numOfDuplication*duplicationCost+ numOfAdditiveTransfers*additiveTransfer+ numOfReplacingTransfers*replacementTransfer;
				
				PrintWriter writer = new PrintWriter(new File(outputFolderPath+"species"+i+"gene"+j+".newick"));
				writer.println(reconciliationCost);
				writer.println(numOfReplacingTransfers);
				writer.close();
				infoScanner.close();
			
			}
		}
	}
	public static RootedTree getGeneTreeFromSimulatedData(String treeString) {
		//System.out.println(treeString.length());
		
	//	System.out.println(treeString);
		boolean flag = false;
		  int length = treeString.length();
		  Stack<Node>stack = new Stack<Node>();
		  for(int i = 0; i < length; i++) {
			  char presentCharacter = treeString.charAt(i);
			  if(presentCharacter=='(' || presentCharacter==',' ) {
				  flag = false;
				  continue;
			  }
			  else if(presentCharacter==')') {
				  Node rightNode = stack.pop();
				  Node leftNode = stack.pop();
				  Node parent = new Node(leftNode.getName()+"$"+rightNode.getName());
				  parent.setLeft(leftNode);
				  parent.setRight(rightNode);
				
				  stack.push(parent);
				  flag = true;
			  }
			  else {
				String name = "";
					name=treeString.substring(i, treeString.indexOf("]",i)+1);
					
				i = treeString.indexOf("]",i)+1;
				if(treeString.charAt(i)==')')i--;
				if(name.contains(":")) {
					String[] partsOfGeneName = name.split(":");
					String geneName = partsOfGeneName[0];
					String secondPart = partsOfGeneName[1];
						
					secondPart = secondPart.substring(secondPart.indexOf('[')+1, secondPart.length()-1);
					String [] descriptiveParts = secondPart.split(" ");
					String id = descriptiveParts[0].split("=")[1];
					String vertexType = descriptiveParts[3].split("=")[1]; 
					
					if(vertexType.equals("Leaf"))
					{	
						geneName = "H"+geneName.split("_")[1];
						Node node = new Node(geneName);
						node.setID(id);
					
						stack.push(node);
					}
				}
				
				
				else {
				Node node = new Node(name);
				if(!flag)
				stack.push(node);
				else flag = false;
				}
			  }
		  }
		  if(stack.isEmpty())
		{
			  System.out.println("Tree: "+treeString.length());
			  System.out.println("Tree completed");
		}
		  Node root = stack.pop();
		  RootedTree geneTree = new RootedTree(root);
		  geneTree.computeAllPreOrderNumbers();
		// FileUtility.setPreOrderNames( geneTree);			
	//	String rangerFormat =  treeToNewick(root);
		
	//	geneTree.setNewickOutputFormat(rangerFormat);
	//	System.out.println(rangerFormat);
	//	System.out.println(geneTree.getLeafSet().size());
		return geneTree;
	}
	private static RootedTree getSpeciesTreeFromSimulatedData(String treeString) {
		//System.out.println(treeString.length());
		
	//	System.out.println(treeString);
		boolean flag = false;
		  int length = treeString.length();
		  Stack<Node>stack = new Stack<Node>();
		  for(int i = 0; i < length; i++) {
			  char presentCharacter = treeString.charAt(i);
			  if(presentCharacter=='(' || presentCharacter==',' ) {
				  flag = false;
				  continue;
			  }
			  else if(presentCharacter==')') {
				  Node rightNode = stack.pop();
				  Node leftNode = stack.pop();
				  Node parent = new Node(leftNode.getName()+"$"+rightNode.getName());
				  parent.setLeft(leftNode);
				  parent.setRight(rightNode);
				
				  stack.push(parent);
				  flag = true;
			  }
			  else {
				String name = "";
					name=treeString.substring(i, treeString.indexOf("]",i)+1);
					
				i = treeString.indexOf("]",i)+1;
				if(treeString.charAt(i)==')')i--;
				if(name.contains(":")) {
					String[] partsOfSpeciesName = name.split(":");
					String speciesName = partsOfSpeciesName[0];
					String secondPart = partsOfSpeciesName[1];
					secondPart = secondPart.substring(secondPart.indexOf('[')+1, secondPart.length()-1);
					String [] descriptiveParts = secondPart.split(" ");
					String id = descriptiveParts[0].split("=")[1];
					String vertexType = descriptiveParts[1].split("=")[1]; 
					
					if(vertexType.equals("Leaf"))
					{	
						Node node = new Node(speciesName);
						node.setID(id);
					
						stack.push(node);
					}
				}
				
				
				else {
				Node node = new Node(name);
				if(!flag)
				stack.push(node);
				else flag = false;
				}
			  }
		  }
		  if(stack.isEmpty())
		{
			  System.out.println("Tree: "+treeString.length());
			  System.out.println("Tree completed");
		}
		  Node root = stack.pop();
		  RootedTree speciesTree = new RootedTree(root);
		 FileUtility.setPreOrderNames( speciesTree);			
	//	String rangerFormat = treeToNewick(speciesTree.getRoot());
	//	speciesTree.setNewickOutputFormat(rangerFormat);
//		System.out.println(rangerFormat);
//		System.out.println(speciesTree.getLeafSet().size());
		return speciesTree;
	}
	public static String  treeToNewick(Node node) {
		
		String output = "";
		String leftPart = "";
		String rightPart = "";
		if(node.getLeft().isLeaf()) {
			leftPart =  node.getLeft().getName();
		}
		else {
			leftPart = treeToNewick(node.getLeft());
		}
		if(node.getRight().isLeaf()) {
			rightPart =  node.getRight().getName();
		}
		else {
			if(node.getLeft().isLeaf())
				rightPart = treeToNewick(node.getRight());
			else
				rightPart = treeToNewick(node.getRight());
				
		}
		output = "("+ leftPart+","+ rightPart+")";
		if(node.isRoot()) {
			
			output+=";";
		}

		return output;
		
	}
}
