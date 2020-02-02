#ifndef VANTAGE_POINT_TREE
#define VANTAGE_POINT_TREE
#include <vector>
#include <iostream>
#include <algorithm>
#include <stack>
#include <tuple>
#include <cmath>
template <typename T>
class vp_tree
{
	class Node
	{
	public:
		Node * inside;
		Node * outside;
		double mu;
		size_t vp;

		Node(Node * in, Node * out, double m, size_t i)
			: inside(in),
			outside(out),
			mu(m),
			vp(i) {}

		~Node() {}
	};

	enum class Node_type {
		left_son, right_son, root
	};

	using point_vec = std::vector<T>;
	using data_vector = std::vector<point_vec>;
	using node_stack = std::stack < std::tuple < Node*, Node_type, bool>>;

	Node * root = nullptr;
	data_vector data;
	size_t current = 0;

	struct Boundaries
	{
		size_t lower;
		size_t upper;
		size_t parent;
		Node_type node_type;
		Boundaries(size_t l, size_t u, size_t p, Node_type le) : lower(l), upper(u), parent(p), node_type(le) {}
	};

	void alloc()
	{
		root = (Node *)operator new (data.size() * sizeof(Node));
	}

	void dealloc()
	{
		operator delete (root);
	}

	void print_vector(const std::vector<T> vec) const
	{
		std::cout << "(";
		for (auto && x : vec)
		{
			std::cout << x << " ";
		}
		std::cout << ")";
	}

		// count appropriate radius to vantage point
	double get_mu(const Boundaries & boundaries, size_t vp) const
	{
		std::vector<double> tmp;
		for (size_t i = boundaries.lower; i < boundaries.upper; i++)
		{
			if (vp != i)
			{
				tmp.push_back(count_distance(data[i], data[vp]));
			}
		}
		return (tmp.size() == 0) ? 0 : find_median(move(tmp));
	}

	double find_median(std::vector<double> && pool) const
	{
		std::nth_element(pool.begin(), pool.begin() + pool.size() / 2, pool.end());
		return pool[pool.size() / 2];
	}

		// Method for updating the buffer of results while searching for k nearest neighbours
	void update_results(std::vector<size_t> & results, double & tau, Node* vp_to_add, const point_vec & point, size_t k) const
	{
		if (results.size() < k)
		{
			results.push_back(vp_to_add->vp);
			if (count_distance(data[vp_to_add->vp], point) > tau)
				tau = count_distance(data[vp_to_add->vp], point);
		}
		else
		{
			if (count_distance(data[vp_to_add->vp], point) > tau)
				return;
			if (k == 1)
			{
				results[0] = vp_to_add->vp;
				tau = count_distance(data[vp_to_add->vp], point);
				return;
			}
			double prev_tau = tau;
			tau = 0;
			bool first = true;
			double dist;
			for (size_t i = 0; i < results.size(); ++i)
			{
				dist = count_distance(data[results[i]], point);
				if (first && prev_tau == dist)
				{
					first = false;
					results[i] = vp_to_add->vp;
				}
				else if (tau < dist)
					tau = dist;
			}
		}
	}
		
		// Method for updating the result buffer during the radius search
	void update_results_radius_search(std::vector<size_t> & results, double radius, Node* vp_to_add, const point_vec & point) const
	{
		if (count_distance(data[vp_to_add->vp], point) > radius)
			return;
		else
			results.push_back(vp_to_add->vp);
	}

	struct current_node
	{
		Node* ptr;
		size_t visited;
		current_node(Node* p, size_t v) : ptr(p), visited(v) {}
	};

		// Search for k nearest neighbours of point, or within the radius around the point
		// I Keep the buffer of results the size of k and value tau which is the distance between the furthest point in buffer and Point
		// If the point is within the `circle`, go to left child (inside) -> same for child (recursion)
		// If tau point is within the `circle` and the left child has already been searched and tau is greater than the distance to the edge of the `circle` around vp,
		// then it is necessary to search right child as well ( same if the point in not within the `circle`)
	std::vector<size_t> search(const point_vec & point, size_t k, double radius) const
	{
		std::vector<size_t> results;
		std::stack< current_node> stack;
		std::stack<Node*> go_up; // to enable forced step to parent
		std::stack< Node*> parents; // pointers to parent nodes
		double tau;
				// always the distance between farthest vp in result tybuffer from point (in the beginning, it is the distance between root & point)
		if (k == 0)
			tau = radius;
		else
			tau = count_distance(data[root->vp], point);
		stack.push(current_node(root, 0));
		parents.push(nullptr);
		while (!stack.empty())
		{
			current_node node = stack.top();
			stack.pop();
			if (!node.ptr)
				return results;
				// update results iff this the first visit of the node
			if (node.visited == 0)
			{
				if (k == 0)
					update_results_radius_search(results, radius, node.ptr, point);
				else
					update_results(results, tau, node.ptr, point, k);
			}
			if (count_distance(data[node.ptr->vp], point) < node.ptr->mu)
				search_inside(go_up, parents, stack, node, tau, point);
			else
				search_outside(go_up, parents, stack, node, tau, point);
		}
		return results;
	}

	void search_inside(std::stack<Node*>& go_up, std::stack<Node*>& parents, std::stack<current_node>& stack, current_node node, double tau, const std::vector<T> & point) const
	{
		if (!go_up.empty() && go_up.top() == node.ptr)
		{		// already searched the children - usually inner nodes, moving to parent 
			go_up.pop();
			Node* parent = parents.top();
			parents.pop();
			stack.push(current_node(parent, 1));
		}	
		else if (node.visited == 0 && node.ptr->inside)
		{		// the first visit of node, going to left child (the point is within mu)
			stack.push(current_node(node.ptr->inside, 0));
			parents.push(node.ptr);
		}
		else if (node.ptr->outside && (node.visited == 1 || !node.ptr->inside) && tau >= node.ptr->mu - count_distance(data[node.ptr->vp], point))
		{		// the latter visit of node, going to right child if results can be found outside mu (condition above)
			stack.push(current_node(node.ptr->outside, 0));
			go_up.push(node.ptr);
			parents.push(node.ptr);
		}
		else //already searched the children - usually leaves or nodes with only one child, moving to parent
		{
			Node* parent = parents.top();
			parents.pop();
			stack.push(current_node(parent, 1));
		}
	}
	 
	void search_outside(std::stack<Node*>& go_up, std::stack<Node*>& parents, std::stack<current_node>& stack, current_node node, double tau, const std::vector<T> & point) const
	{
			//already searched the children - usually inner nodes, moving to parent 
		if (!go_up.empty() && go_up.top() == node.ptr)
		{
			go_up.pop();
			Node* parent = parents.top();
			parents.pop();
			stack.push(current_node(parent, 1));
		}
		else if (node.visited == 0 && node.ptr->outside)
		{		// the first visit of node, going to right child (the point is not within mu)
			stack.push(current_node(node.ptr->outside, 0));
			parents.push(node.ptr);
		}
		else if (node.ptr->inside && (node.visited == 1 || !node.ptr->outside) && tau >= count_distance(data[node.ptr->vp], point) - node.ptr->mu)
		{		// the latter visit of node, going to left child if results can be found within mu (condition above)
			stack.push(current_node(node.ptr->inside, 0));
			go_up.push(node.ptr);
			parents.push(node.ptr);
		}
		else //already searched the children - usually leaves or nodes with only one child, moving to parent
		{
			Node* parent = parents.top();
			parents.pop();
			stack.push(current_node(parent, 1));
		}
	}

	size_t choose_vp(const Boundaries & boundaries) const
	{
		return (((size_t)rand() %( boundaries.upper - boundaries.lower)) +boundaries.lower);
	}

	void swap(vp_tree & other) 
	{
		std::swap(data, other.data);
		std::swap(current, other.current);
		std::swap(root, other.root);
	}

		// constructs the whole vantage point tree from data 
	void build_tree()
	{
		std::stack<Boundaries> stack;
		size_t current = 0;
		stack.emplace(0, data.size(), 0, Node_type::root);

		while (!stack.empty())
		{
			Boundaries boundaries = stack.top();
			stack.pop();
			size_t vp = choose_vp(boundaries);
			double mu = get_mu(boundaries, vp);
			size_t half = boundaries.lower + (boundaries.upper - boundaries.lower) / 2;
			std::swap(data[vp], data[half]);
			vp = half;
			correct_misplaced_points(boundaries.lower, half + 1, vp, mu, boundaries);
			insert_node(mu, current, vp, boundaries, stack);
			++current;
		}
	}

	void insert_node(double mu, size_t current, size_t vp, const Boundaries & boundaries, std::stack<Boundaries> & stack)
	{
		new (root + current) Node(nullptr, nullptr, mu, vp);
			// left (inside) son
		if (vp != 0 && vp > boundaries.lower)
			stack.emplace(boundaries.lower, vp, current, Node_type::left_son);
			// right (outside) son
		if (vp + 1 < boundaries.upper)
			stack.emplace(vp + 1, boundaries.upper, current, Node_type::right_son);
		
			// Add the child node(this node) to parent node created earlier
		if (boundaries.node_type == Node_type::left_son)
			(root + boundaries.parent)->inside = (root + current);
		else if (boundaries.node_type == Node_type::right_son)
			(root + boundaries.parent)->outside = (root + current);
	}

	void correct_misplaced_points(size_t in, size_t out, size_t vp, double mu, const Boundaries & boundaries)
	{
			// I dont know beforehand if the points on the edge of `circle` of vp should be inside/outside
			// therefore 2x while
		size_t _in = in, _out = out;
		while (in < vp && out < boundaries.upper)
		{		// the outside points are correctly arranged
			while (in < vp && count_distance(data[in], data[vp]) < mu)
				++in;
			while (out < boundaries.upper && count_distance(data[out], data[vp]) >= mu)
				++out;
			if (in < vp && out < boundaries.upper)
				std::swap(data[in], data[out]);
			if (in < vp)
				++in;
			if (out < boundaries.upper)
				++out;
		}
		in = _in;
		out = _out;
			
		while (in < vp && out < boundaries.upper)
		{		// the inside points are correctly arranged
			while (in < vp && count_distance(data[in], data[vp]) <= mu)
				++in;
			while (out < boundaries.upper && count_distance(data[out], data[vp]) > mu)
				++out;
			if (in < vp && out < boundaries.upper)
				std::swap(data[in], data[out]);
			if (in < vp)
				++in;
			if (out < boundaries.upper)
				++out;
		}

	}

public:
	vp_tree(const data_vector & other)
	{
		data = other;
		alloc();
		build_tree();
	}

	vp_tree(data_vector && other)
	{
		data = move(other);
		alloc();
		build_tree();
	}

	vp_tree(vp_tree && other)
	{
		root = std::move(other.root);
		data = std::move(other.data);
		current = std::move(current);
		other.root = nullptr;
	}

	vp_tree(const vp_tree & other)
	{
			// need deep copy 
		data.resize(other.data.size());
		for (size_t i = 0; i < other.data.size(); ++i)
		{
			point_vec vec{};
			vec.resize(other.data[i].size());
			for (size_t j = 0; j < other.data[i].size(); ++j)
			{
				vec[j] = other.data[i][j];
			}
			data[i] = std::move(vec);
		}
		alloc();
		build_tree();
	}

	vp_tree & operator = (const vp_tree & other)
	{
		vp_tree tmp(other);
		swap(tmp);
		return *this;
	}

	vp_tree & operator = (vp_tree && other)
	{
		dealloc();
		root = std::move(other.root);
		data = std::move(other.data);
		current = std::move(current);
		other.root = nullptr;
		return *this;
	}

	~vp_tree()
	{
		if (root)
			dealloc();
	}

		// if k < #of point in vp, then UB
	data_vector search_k_nearest_neighbours(const point_vec & point, size_t k) const
	{
		if (k == 0)
			return data_vector();
		std::vector<size_t> res = search(point, k, 0);
		data_vector vecs;
		for (size_t i = 0; i < res.size(); i++)
		{
			vecs.push_back(data[res[i]]);
		}
		return vecs;
	}

	data_vector search_within_radius(const point_vec & point, double radius) const
	{
		std::vector<size_t> res = search(point, 0, radius);
		data_vector vecs;
		for (size_t i = 0; i < res.size(); i++)
		{
			vecs.push_back(data[res[i]]);
		}
		return vecs;
	}

		// Override to change the metric
	double virtual count_distance(const point_vec & a, const point_vec & b) const
	{
		double diff = 0;
		for (size_t i = 0; i < a.size(); ++i)
		{
			diff += (a[i] - b[i]) * (a[i] - b[i]);
		}
		return std::sqrt(diff);
	}

	void print_data()
	{
		for (size_t i = 0; i < data.size(); ++i)
		{
			print_vector(data[i]);
			std::cout << std::endl;
		}
	}

	void print_tree()
	{
		for (size_t i = 0; i < data.size(); i++)
		{
			std::cout << std::endl;
			if ((root + i)->inside != nullptr)
				std::cout << "inside(left_son)  : " << (root + i)->inside << " >>>>> (" << data[static_cast<Node*>((root + i)->inside)->vp][0] << ',' << data[static_cast<Node*>((root + i)->inside)->vp][1] << ')' << std::endl;
			if ((root + i)->outside != nullptr)
				std::cout << "outside(right_son): " << (root + i)->outside << " >>>>> (" << data[static_cast<Node*>((root + i)->outside)->vp][0] << ',' << data[static_cast<Node*>((root + i)->outside)->vp][1] << ')' << std::endl;
			std::cout << "mu(radius)        : " << (root + i)->mu << std::endl;
			std::cout << "vantage point     : " << (root + i)->vp << " >>>>> (" << data[(root + i)->vp][0] << ',' << data[(root + i)->vp][1] << ')' << std::endl;
			std::cout << "***************************************";
		}
	}
};
#endif
