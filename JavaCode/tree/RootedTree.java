package tree;

import java.io.IOException;
import java.util.*;

import edge.Edge;
import node.Node;
import utility.FileUtility;

public class RootedTree {
	Node root;
	String type;
	boolean dated;
	String newickOutputFormat;
	ArrayList<Node> allNodesByPreOrder = new ArrayList<Node>();
	ArrayList<Node> allNodesByPostOrder = new ArrayList<Node>();
	ArrayList<Node> allNodesByInOrder = new ArrayList<Node>();
	ArrayList<Edge> allBranches = new ArrayList<Edge>();
	ArrayList<Node> leafSet = new ArrayList<Node>();
	ArrayList<Node> allPostOrderInternalNodes = new ArrayList<Node>();
	ArrayList<Node> allPreOrderInternalNodes = new ArrayList<Node>();
	
	public RootedTree(Node root) {
		this.root = root;
	}

	public RootedTree() {
	}

	public Node getRoot() {
		return root;
	}

	public String getNewickOutputFormat() {
		return newickOutputFormat;
	}

	public void setNewickOutputFormat(String newickOutputFormat) {
		this.newickOutputFormat = newickOutputFormat;
	}

	
	public ArrayList<Node> getLeafSet() {
		return leafSet;
	}

	public void setLeafSet(ArrayList<Node> leafSet) {
		this.leafSet = leafSet;
	}

	public ArrayList<Node> getAllPostOrderInternalNodes() {
		return allPostOrderInternalNodes;
	}

	public void setAllPostOrderInternalNodes(ArrayList<Node> allPostOrderInternalNodes) {
		this.allPostOrderInternalNodes = allPostOrderInternalNodes;
	}

	public ArrayList<Node> getAllPreOrderInternalNodes() {
		return allPreOrderInternalNodes;
	}

	public void setAllPreOrderInternalNodes(ArrayList<Node> allPreOrderInternalNodes) {
		this.allPreOrderInternalNodes = allPreOrderInternalNodes;
	}

	public ArrayList<Edge> getAllBranches() {
		return allBranches;
	}

	public void setAllBranches(ArrayList<Edge> allBranches) {
		this.allBranches = allBranches;
	}

	public void setRoot(Node root) {
		this.root = root;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public boolean isDated() {
		return dated;
	}

	public void setDated(boolean dated) {
		this.dated = dated;
	}

	public ArrayList<Node> getAllNodesByPreOrder() {
		return allNodesByPreOrder;
	}

	public void setAllNodesByPreOrder(ArrayList<Node> allNodesByPreOrder) {
		this.allNodesByPreOrder = allNodesByPreOrder;
	}

	public ArrayList<Node> getAllNodesByPostOrder() {
		return allNodesByPostOrder;
	}

	public void setAllNodesByPostOrder(ArrayList<Node> allNodesByPostOrder) {
		this.allNodesByPostOrder = allNodesByPostOrder;
	}

	public ArrayList<Node> getAllNodesByInOrder() {
		return allNodesByInOrder;
	}

	public void setAllNodesByInOrder(ArrayList<Node> allNodesByInOrder) {
		this.allNodesByInOrder = allNodesByInOrder;
	}

	public static void main(String[] args) throws IOException {
		Node a = new Node("a");
		Node b = new Node("b");
		Node c = new Node("c");
		Node d = new Node("d");
		Node e = new Node("e");
		Node x = new Node("x");
		Node y = new Node("y");
		Node z = new Node("z");
		Node w = new Node("w");

		w.setLeft(y);
		w.setRight(z);

		y.setLeft(x);
		y.setRight(c);

		z.setLeft(d);
		z.setRight(e);

		x.setLeft(a);
		x.setRight(b);
		RootedTree tree = new RootedTree(w);
		tree.computeAllInOrderNumbers();
		tree.computeAllPreOrderNumbers();
		tree.computeAllPostOrderNumbers();

		System.out.println(a);
		System.out.println(b);
		System.out.println(c);
		System.out.println(d);
		System.out.println(e);
		System.out.println(x);
		System.out.println(y);
		System.out.println(z);
		System.out.println(w);

		System.out.println("Left most of y: " + y.getCalculatedLeftMostDescendent());
		System.out.println("Right most of w: " + w.getCalculatedRightMostDescendent());
		System.out.println(FileUtility.getSpeciesTreeFromFile("input.txt"));
	}

	public ArrayList<Node> computeLeafSet() {
		//ArrayList<Node> list = new ArrayList<Node>();
		Stack<Node> stack = new Stack<Node>();
		stack.push(root);

		while (!stack.isEmpty()) {
			Node node = stack.pop();

			if (node.isLeaf())
				leafSet.add(node);

			if (node.getRight() != null) {
				stack.push(node.getRight());
			}
			if (node.getLeft() != null) {
				stack.push(node.getLeft());
			}
		}

		return leafSet;
	}

	public void computeAllPreOrderNumbers() {
	
		Stack<Node> stack = new Stack<Node>();
		stack.push(root);
		root.setDepth(0);
		int number = 0;
		// System.out.print("Pre-Order: ");
		while (!stack.isEmpty()) {
			Node node = stack.pop();
			if (node.getParent() != null) {

				Node parent = node.getParent();
				node.setDepth(parent.getDepth() + 1);
			}
			node.setPreOrderNumber(number);
			allNodesByPreOrder.add(node);
			// System.out.print(node.getName()+" ");
			number++;

			if (node.getRight() != null) {
				stack.push(node.getRight());
			}
			if (node.getLeft() != null) {
				stack.push(node.getLeft());
			}
		}
		// System.out.println();
	}

	public void computeAllPostOrderNumbers() {
		Stack<Node> s1 = new Stack<Node>();
		Stack<Node> s2 = new Stack<Node>();

		if (root != null)
			s1.push(root);

		// Run while first stack is not empty
		while (!s1.isEmpty()) {
			// Pop an item from s1 and push it to s2
			Node temp = s1.pop();
			s2.push(temp);

			// Push left and right children of
			// removed item to s1
			if (temp.getLeft() != null)
				s1.push(temp.getLeft());
			if (temp.getRight() != null)
				s1.push(temp.getRight());
		}
		int number = 0;
		// Print all elements of second stack
		// System.out.print("Post-Order: ");
		while (!s2.isEmpty()) {
			Node temp = s2.pop();
			allNodesByPostOrder.add(temp);
			// System.out.print(temp.getName()+" ");
			temp.setPostOrderNumber(number);
			number++;
		}
		// System.out.println();
	}

	public void computeAllInOrderNumbers() {
		Stack<Node> s = new Stack<Node>();
		Node current = root;
		int number = 0;
		// System.out.print("In-Order: ");
		// traverse the tree
		while (current != null || s.size() > 0) {

			/*
			 * Reach the left most Node of the curr Node
			 */
			while (current != null) {
				/*
				 * place pointer to a tree node on the stack before traversing the node's left
				 * subtree
				 */
				s.push(current);
				current = current.getLeft();
			}

			/* Current must be NULL at this point */
			current = s.pop();
			allNodesByInOrder.add(current);
			// System.out.print(current.getName()+" ");
			current.setInOrderNumber(number);
			number++;
			/*
			 * we have visited the node and its left subtree. Now, it's right subtree's turn
			 */
			current = current.getRight();
		}
		// System.out.println();
	}

	public ArrayList<Node> computeAllPostOrderInternalNodes() {
		int size = allNodesByPostOrder.size();
		for (int i = 0; i < size; i++) {
			Node node = allNodesByPostOrder.get(i);
			if (!node.isLeaf())
				allPostOrderInternalNodes.add(node);
		}
		return allPostOrderInternalNodes;
	}

	public ArrayList<Node> computeAllPreOrderInternalNodes() {
		int size = allNodesByPreOrder.size();
		for (int i = 0; i < size; i++) {
			Node node = allNodesByPreOrder.get(i);
			if (!node.isLeaf())
				allPreOrderInternalNodes.add(node);
		}
		return allPreOrderInternalNodes;
	}

	public void setAllBranches() {
		ArrayList<Node> allInternalNodes = allPreOrderInternalNodes;
		for (Node node : allInternalNodes) {
			node.setBranches();
			allBranches.add(node.getLeftBranch());
			allBranches.add(node.getRightBranch());
		}
	}
	
}
