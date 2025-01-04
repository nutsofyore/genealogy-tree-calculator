/*
Authors: Yaprak Deniz Yurt, Virani Baskurt

The code is generated for my doctoral dissertation. It aimed to calculate the novelty and variety scores of sketches automatically and error-free. 
The keywords were generated from the sketches of subjects following the FBS (function-behavior-structure) ontology and then transferred into the code space.
The code calculates these two metrics based on the formulas for ideation effectiveness: Hernandez, Shah, and Smith (2010), Nelson et al.(2009), and Bayirli (2018). 
After reviewing the outputs, I have continued with Bayirli's approach.

The resources for the formulas:
Bayirli, U. (2018). FICTIONATION IDEA GENERATION TOOL FOR PRODUCT DESIGN EDUCATION UTILIZING WHAT-IF SCENARIOS OF DESIGN FICTION: A MIXED METHOD STUDY A THESIS SUBMITTED TO THE GRADUATE SCHOOL OF NATURAL AND APPLIED SCIENCES OF MIDDLE EAST TECHNICAL UNIVERSITY.
Hernandez, N. V., Shah, J. J., & Smith, S. M. (2010). Understanding design ideation mechanisms through multilevel aligned empirical studies. Design Studies, 31(4), 382-410. 
https://doi.org/10.1016/j.destud.2010.04.001
Nelson, B. A., Wilson, J. O., Rosen, D., & Yen, J. (2009). Refined metrics for measuring ideation effectiveness. Design Studies, 30(6), 737-743. 
https://doi.org/10.1016/j.destud.2009.07.002

Created for: Yurt, Y. D. (2025). Investigating the Role of Incubation in Design Creativity. [Unpublished doctoral dissertation]. Middle East Technical University.
*/

#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <queue>
#include <stack>
#include <functional>

const std::string P1 = "P1";
const std::string P2 = "P2";
const std::string P3 = "P3";
const std::string P4 = "P4";
const std::string P5 = "P5";
const std::string P6 = "P6";
const std::string P7 = "P7";
const std::string P8 = "P8";
const std::string P9 = "P9";
const std::string P10 = "P10";
const std::string P11 = "P11";
const std::string P12 = "P12";
const std::string P13 = "P13";
const std::string P14 = "P14";
const std::string P15 = "P15";
const std::string P16 = "P16";
const std::string P17 = "P17";
const std::string P18 = "P18";
const std::string P19 = "P19";
const std::string P20 = "P20";
const std::string P21 = "P21";

// List of participants
const std::vector<std::string> participants{
    P1,
    P2,
    P3,
    P4,
    P5,
    P6,
    P7,
    P8,
    P9,
    P10,
    P11,
    P12,
    P13,
    P14,
    P15,
    P16,
    P17,
    P18,
    P19,
    P20,
    P21,
};

struct Node
{
    std::string name;
    int level = -1; // Level of the node in the hierarchy
    std::vector<std::shared_ptr<Node>> children; // Child nodes

    Node(std::string name)
        : name(name) {}

    Node(std::string name, int level)
        : name(name), level(level) {}
};

class Tree
{
public:
    std::shared_ptr<Node> root;
    int depth;
    Tree(const std::string rootName) : root(std::make_shared<Node>(rootName)), depth(0)
    {
        initializeLevels(this->root);
    }

    // Adds a new node to the tree under a specified parent
    bool addNode(const std::string &nodeName, const std::string &parentName)
    {
        if (getNode(nodeName))
        {
            std::cout << nodeName << " ALREADY EXIST! FAILED TO ADD" << std::endl;
            return false;
        }

        std::shared_ptr<Node> parentNode = getNode(parentName);
        if (!parentNode)
        {
            std::cout << parentName << " PARENT DOESN'T EXIST! FAILED TO ADD" << std::endl;

            return false;
        }

        std::shared_ptr<Node> addingNode = std::make_shared<Node>(nodeName);
        parentNode->children.emplace_back(std::move(addingNode));
        initializeLevels(root);
        return true;
    }

    std::vector<int> getNumberOfNodesInEachLevel()
    {
        std::vector<int> result(depth);
        auto getLevel = [&result](const std::shared_ptr<Node> &node)
        {
            result[node->level]++;
        };

        bfs(root, getLevel);

        return result;
    }

    std::shared_ptr<Tree> getCopyTreeForNodes(const std::vector<std::string> &searchingNodes)
    {
        std::shared_ptr<Node> copyRoot = getTreeForNodesInternal(root, searchingNodes);
        if (!copyRoot)
        {
            return NULL;
        }
        // TODO this method copies the tree by duplicating
        std::shared_ptr<Tree> copyTree = std::make_shared<Tree>(copyRoot->name);

        auto getLevel = [&copyTree](const std::shared_ptr<Node> &node)
        {
            for (const auto &child : node->children)
            {
                copyTree->addNode(child->name, node->name);
            }
        };

        bfs(copyRoot, getLevel);

        return copyTree;
    }

    // source: https://stackoverflow.com/a/51730733/5225074
    void print()
    {
        std::function<void(const std::shared_ptr<Node> &, const std::string &, bool)> print_internal =
            [&print_internal](const std::shared_ptr<Node> &node, const std::string &prefix, bool isFirstChild) -> void
        {
            if (node != nullptr)
            {
                std::cout << prefix;
                std::cout << (isFirstChild ? "|__ " : "|__ ");
                std::cout << node->name << std::endl;
                for (size_t i = 0; i < node->children.size(); i++)
                {
                    print_internal(node->children[i], prefix + (isFirstChild ? "|\t" : "\t"), i == 0);
                }
            }
        };

        print_internal(root, "", false);
    }

private:
    // Performs breadth-first search (BFS) to apply a function to each node
    void bfs(const std::shared_ptr<Node> &root,
             std::function<void(const std::shared_ptr<Node> &)> predicate)
    {
        std::queue<std::shared_ptr<Node>> searchQueue;
        searchQueue.emplace(root);

        while (!searchQueue.empty())
        {
            auto currentNode = searchQueue.front();
            searchQueue.pop();

            predicate(currentNode);

            for (const auto &child : currentNode->children)
            {
                searchQueue.emplace(child);
            }
        }
    }
    
    // Performs depth-first search (DFS) to apply a function to each node
    void dfs(const std::shared_ptr<Node> &root,
             std::function<bool(const std::shared_ptr<Node> &)> predicate = nullptr)
    {

        std::stack<std::shared_ptr<Node>> searchStack;
        searchStack.emplace(root);
        while (!searchStack.empty())
        {
            auto currentNode = searchStack.top();
            searchStack.pop();

            predicate(currentNode);

            for (const auto &child : currentNode->children)
            {
                searchStack.emplace(child);
            }
        }
    }

    std::shared_ptr<Node> getTreeForNodesInternal(const std::shared_ptr<Node> &node, const std::vector<std::string> &searchingNodes)
    {
        // TODO this function has a bug. If you're looking for both c1 and l1, where c1 is the parent for l1, it only returns c1. It doesn't go any further to find l1
        if (std::find(searchingNodes.begin(), searchingNodes.end(), node->name) != searchingNodes.end())
        {
            return std::make_shared<Node>(node->name);
        }

        std::shared_ptr<Node> copyNode = nullptr;
        for (const auto &child : node->children)
        {
            auto copyChild = getTreeForNodesInternal(child, searchingNodes);
            if (copyChild)
            {
                if (!copyNode)
                    copyNode = std::make_shared<Node>(node->name, node->level);
                copyNode->children.emplace_back(copyChild);
            }
        }

        return copyNode;
    }

    // Retrieves a specific node by name
    std::shared_ptr<Node> getNode(const std::string &nodeName)
    {
        std::shared_ptr<Node> foundNode = nullptr;
        auto getByName = [&nodeName, &foundNode](const std::shared_ptr<Node> &node)
        {
            if (node->name == nodeName)
            {
                foundNode = node;
            }
        };

        bfs(root, getByName);
        return foundNode;
    }

    // Assigns level values to each node in the tree
    void initializeLevels(const std::shared_ptr<Node> &node)
    {
        int currentLevel = 0;
        int treeDepth = -1;

        std::function<void(const std::shared_ptr<Node> &, int)> initializeLevels_internal =
            [&initializeLevels_internal, &treeDepth](const std::shared_ptr<Node> &node, int currentLevel)
        {
            node->level = currentLevel;
            ++currentLevel;
            treeDepth = std::max(treeDepth, currentLevel);
            for (const auto &child : node->children)
            {
                initializeLevels_internal(child, currentLevel);
            }
        };

        initializeLevels_internal(root, 0);
        depth = treeDepth;
    }
};

std::vector<std::string> getLeafNamesForParticipant(const std::string &participant,
                                                    const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    std::vector<std::string> result;

    for (auto it = leafParticipantMap.cbegin(); it != leafParticipantMap.cend(); ++it)
    {
        if (std::find(it->second.cbegin(), it->second.cend(), participant) != it->second.cend() &&
            std::find(result.cbegin(), result.cend(), it->first) == result.cend())
        {
            result.push_back(it->first);
        }
    }
    if (result.size() == 0)
    {
        std::cout << participant << " DOESN'T EXIST IN THE MAP" << std::endl;
    }

    return result;
}

std::shared_ptr<Tree> getTreeForParticipant(const std::shared_ptr<Tree> &tree, const std::string &participant,
                                            const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    const std::vector<std::string> searchingNodes = getLeafNamesForParticipant(participant, leafParticipantMap);
    return tree->getCopyTreeForNodes(searchingNodes);
}

enum class VarietyScoreCalculationMode
{
    Nelson /*substracts 1*/,
    Shah
};

float calculateTotalNumberOfIdeas(const std::shared_ptr<Node> root,
                                  const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    if (root->children.size() == 0)
    {
        auto result = leafParticipantMap.find(root->name);
        if (result == leafParticipantMap.end())
        {
            std::cout << root->name << " DOESN'T EXIST IN THE LEAF PARTICIPANT MAP" << std::endl;
        }
        return result->second.size();
    }

    float numberOfIdeas = 0;
    for (const auto &child : root->children)
    {
        numberOfIdeas += calculateTotalNumberOfIdeas(child, leafParticipantMap);
    }

    return numberOfIdeas;
}

float calculateVarietyScoreShah(
    const std::shared_ptr<Tree> &tree,
    const std::vector<int> &levelWeights,
    const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    float totalNumberOfIdeas = calculateTotalNumberOfIdeas(tree->root, leafParticipantMap);

    float score = 0;
    auto numberOfNodesInEachLevel = tree->getNumberOfNodesInEachLevel();
    for (int level = 0; level < numberOfNodesInEachLevel.size(); level++)
    {
        int numberOfNode = numberOfNodesInEachLevel[level];
        score += numberOfNode * levelWeights[level];
    }

    float result = score / totalNumberOfIdeas;
    return result;
}

float calculateVarietyScoreNelson(
    const std::shared_ptr<Tree> &tree,
    const std::vector<int> &levelWeights,
    const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    float totalNumberOfIdeas = calculateTotalNumberOfIdeas(tree->root, leafParticipantMap) - 1.0f;

    float score = 0;
    auto numberOfNodesInEachLevel = tree->getNumberOfNodesInEachLevel();
    for (int level = 0; level < numberOfNodesInEachLevel.size(); level++)
    {
        int numberOfNode = numberOfNodesInEachLevel[level] - 1;
        if (numberOfNode < 0)
        {
            std::cout << "numberOfNode is negative: " << numberOfNode << std::endl;
        }
        score += numberOfNode * levelWeights[level];
    }

    float result = score / totalNumberOfIdeas;
    return result;
}

float calculateVarietyScoreBayirliShah(
    const std::shared_ptr<Tree> &tree,
    const std::vector<int> &levelWeights,
    float treeScore)
{

    float score = 0;
    auto numberOfNodesInEachLevel = tree->getNumberOfNodesInEachLevel();
    for (int level = 0; level < numberOfNodesInEachLevel.size(); level++)
    {
        int numberOfNode = numberOfNodesInEachLevel[level];
        score += numberOfNode * levelWeights[level];
    }

    float result = score * 10.0f / treeScore;
    return result;
}

float calculateVarietyScoreBayirliNelson(
    const std::shared_ptr<Tree> &tree,
    const std::vector<int> &levelWeights,
    float treeScore)
{
    float score = 0;
    auto numberOfNodesInEachLevel = tree->getNumberOfNodesInEachLevel();
    for (int level = 0; level < numberOfNodesInEachLevel.size(); level++)
    {
        int numberOfNode = numberOfNodesInEachLevel[level] - 1;
        if (numberOfNode < 0)
        {
            std::cout << "numberOfNode is negative: " << numberOfNode << std::endl;
        }
        score += numberOfNode * levelWeights[level];
    }

    float result = score * 10.0f / treeScore;
    return result;
}

float calculateTreeScore(const std::shared_ptr<Tree> &tree,
                         const std::vector<int> &levelWeights)
{
    float score = 0;
    auto numberOfNodesInEachLevel = tree->getNumberOfNodesInEachLevel();
    for (int level = 0; level < numberOfNodesInEachLevel.size(); level++)
    {
        int numberOfNode = numberOfNodesInEachLevel[level];
        score += numberOfNode * levelWeights[level];
    }

    return score;
}

float calculateNoveltyScoreForParticipant(const std::string &participant,
                                          const std::unordered_map<std::string, float> &leafNoveltyScoreMap,
                                          const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    float score = 0;
    const auto leaves = getLeafNamesForParticipant(participant, leafParticipantMap);
    for (const auto &leaf : leaves)
    {
        float leafNoveltyScore = leafNoveltyScoreMap.at(leaf);
        // float numberOfIdeasForStudent = 1.0f; //we ignore it
        score += leafNoveltyScore;
    }
    score = score / leaves.size();

    return score;
}

float getNumberOfIdeasInALeaf(const std::string &leafName,
                              const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    auto result = leafParticipantMap.find(leafName);
    if (result == leafParticipantMap.end())
    {
        std::cout << "getNumberOfIdeasInALeaf: " << leafName << " DOESN'T EXIST IN THE LEAF PARTICIPANT MAP" << std::endl;
    }
    return result->second.size();
}

float calculatNoveltyScoreForLeaf(const std::string &leafName,
                                  float treeNumberOfIdeas,
                                  const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    // N=(T-S)*10/T
    // T tree total number of ideas
    // S number of ideas for leaf
    float s = getNumberOfIdeasInALeaf(leafName, leafParticipantMap);
    float t = treeNumberOfIdeas;
    float result = (t - s) * 10.0f / t;
    return result;
}

std::unordered_map<std::string, float> calculateNoveltyScoreForLeaves(float treeTotalNumberOfIdeas,
                                                                      const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    std::unordered_map<std::string, float> leafVarietyScoreMap;
    for (const auto &leafParticipant : leafParticipantMap)
    {
        const auto &leafName = leafParticipant.first;
        leafVarietyScoreMap[leafName] = calculatNoveltyScoreForLeaf(leafName, treeTotalNumberOfIdeas, leafParticipantMap);
    }
    return leafVarietyScoreMap;
}

std::unordered_map<std::string, float> calculateNoveltyScoreForParticipants(const std::unordered_map<std::string, float> &leafNoveltyScoreMap,
                                                                            const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    auto checkLeafParticipantMapContains = [&leafParticipantMap](const std::string &participant) -> bool
    {
        for (const auto &pair : leafParticipantMap)
        {
            if (std::find(pair.second.begin(), pair.second.end(), participant) != pair.second.end())
            {
                return true;
            }
        }
        return false;
    };

    std::unordered_map<std::string, float> participantNoveltyScoreMap;
    for (const auto &participant : participants)
    {
        if (checkLeafParticipantMapContains(participant))
        {
            participantNoveltyScoreMap[participant] = calculateNoveltyScoreForParticipant(participant, leafNoveltyScoreMap, leafParticipantMap);
        }
    }
    return participantNoveltyScoreMap;
}

std::unordered_map<std::string, const std::vector<std::string>> getLeafMapForParticipant(
    const std::string &participant,
    const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    std::unordered_map<std::string, const std::vector<std::string>> filteredMap;

    for (const auto &entry : leafParticipantMap)
    {
        const auto &participants = entry.second;

        if (std::find(participants.begin(), participants.end(), participant) != participants.end())
        {
            filteredMap.emplace(entry.first, std::vector<std::string>{participant});
        }
    }
    return filteredMap;
}

void printNoveltyScores(const std::string &treeName, const std::shared_ptr<Tree> &tree,
                        const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    float totalNumberOfIdeas = calculateTotalNumberOfIdeas(tree->root, leafParticipantMap);
    auto leafNoveltyScoresMap = calculateNoveltyScoreForLeaves(totalNumberOfIdeas, leafParticipantMap);
    auto participantNoveltyScoreMap = calculateNoveltyScoreForParticipants(leafNoveltyScoresMap, leafParticipantMap);

    std::cout << "--" << treeName << " NOVELTY SCORES--" << ::std::endl;
    float sum = 0.0f;
    for (const auto &participantScore : participantNoveltyScoreMap)
    {
        sum += participantScore.second;
        std::cout << participantScore.first << ": " << participantScore.second << std::endl;
    }
    std::cout << treeName << "novelty score average: " << sum / participantNoveltyScoreMap.size() << ::std::endl;
    std::cout << ::std::endl;
}

void printVarietyScores(const std::string &treeName,
                        const std::shared_ptr<Tree> &tree,
                        const std::vector<int> &levelWeights,
                        const std::unordered_map<std::string, const std::vector<std::string>> &leafParticipantMap)
{
    std::cout << "--" << treeName << "VARIETY SCORES--" << ::std::endl;

    float treeScore = calculateTreeScore(tree, levelWeights);
    std::cout << "Tree Score " << treeScore << std::endl;
    for (const auto &participant : participants)
    {
        const auto participantTree = getTreeForParticipant(tree, participant, leafParticipantMap);
        if (participantTree)
        {
            const auto participantLeafMap = getLeafMapForParticipant(participant, leafParticipantMap);

            float pVarietyScoreNelson = calculateVarietyScoreNelson(participantTree, levelWeights, participantLeafMap);
            float pVarietyScoreShah = calculateVarietyScoreShah(participantTree, levelWeights, participantLeafMap);
            float pVarietyScoreBayirliNelson = calculateVarietyScoreBayirliNelson(participantTree, levelWeights, treeScore);
            float pVarietyScoreBayirliShah = calculateVarietyScoreBayirliShah(participantTree, levelWeights, treeScore);

            std::cout << "For " << participant << " Nelson: " << pVarietyScoreNelson << " Shah: " << pVarietyScoreShah << " BayirliNelson: " << pVarietyScoreBayirliNelson << " BayirliShah: " << pVarietyScoreBayirliShah << ::std::endl;
        }
    }
    std::cout << ::std::endl;
}

void printTestResult(const std::string label, float calculated, float expected)
{
    const std::string successful = "\033[1;32mSUCCESSFUL\033[0m";
    const std::string fail = "\033[1;31mFAIL\033[0m";
    std::cout
        << label << ": "
        << "Calculated: " << calculated
        << " Expected: " << expected << " Result: "
        << ((calculated == expected) ? successful : fail)
        << std::endl;
}

void CalculateGenealogyTreeScore()
{
    std::vector<int> levelWeights = {0, 25, 20, 15, 10, 5, 3, 1};

#pragma region NODE_NAMES
    // Make sure the node names are unique
    const std::string Ideas = "Ideas";

    const std::string BodyEnclosure = "BodyEnclosure";

    const std::string Form = "Form";
    const std::string Dimensions = "Dimensions";
    const std::string Colors = "Colors";
    const std::string Orientation = "Orientation";

    const std::string Rigid1 = "Rigid1";
    const std::string Flexible1 = "Flexible1";
    const std::string SelfStanding = "SelfStanding";
    const std::string VerticalOrientation = "VerticalOrientation";
    const std::string HorizontalOrientation = "HorizontalOrientation";
    const std::string Hangable = "Hangable";

    const std::string Prismatic = "Prismatic";
    const std::string Conical = "Conical";
    const std::string SemiSpherical = "SemiSpherical";
    const std::string Cylindrical = "Cylindrical";
    const std::string PurseLike = "PurseLike";
    const std::string Sculpted = "Sculpted";
    const std::string Wrappable = "Wrappable";
    const std::string BeanBag = "BeanBag";

    const std::string Rectangular = "Rectangular";
    const std::string LShaped = "LShaped";
    const std::string Triangular = "Triangular";
    const std::string Hexagon = "Hexagon";
    const std::string Hourglass = "Hourglass";
    const std::string Assymetrical = "Assymetrical";
    const std::string Straight = "Straight";
    const std::string SemiCylindrical = "SemiCylindrical";
    const std::string ForOtherFunctions = "ForOtherFunctions";
    const std::string ToMountAStrap = "ToMountAStrap";
    const std::string ToGraspBody = "ToGraspBody";
    const std::string ToFormFeet = "ToFormFeet";

    const std::string Multisurface = "Multisurface";
    const std::string CurvedSurfaces = "CurvedSurfaces";
    const std::string FlatSurfaces = "FlatSurfaces";

#pragma endregion NODE_NAMES

#pragma region LEAF_PARTICIPANT_MAP
    std::unordered_map<std::string, const std::vector<std::string>> leafParticipantMap =
        {
            {Dimensions, {P16}},
            {Colors, {P10, P20}},

            {SelfStanding, {P4, P6, P8, P10, P12, P14, P16, P18, P20}},
            {VerticalOrientation, {P8, P12, P14, P16, P18, P20}},
            {HorizontalOrientation, {P10, P14, P16, P18}},
            {Hangable, {P10, P12, P16}},

            {Conical, {P20}},
            {SemiSpherical, {P10}},
            {PurseLike, {P10}},
            {Wrappable, {P16}},
            {BeanBag, {P20}},

            {LShaped, {P6}},
            {Triangular, {P14}},
            {Hexagon, {P12, P18}},
            {Hourglass, {P8, P20}},
            {Assymetrical, {P6}},
            {Straight, {P8, P12, P18}},
            {SemiCylindrical, {P16}},
            {ForOtherFunctions, {P8}},
            {ToMountAStrap, {P6}},
            {ToGraspBody, {P12, P20}},
            {ToFormFeet, {P6, P10}},

            {Multisurface, {P10}},
            {CurvedSurfaces, {P4, P12, P14}},
            {FlatSurfaces, {P4, P8}},
        };
#pragma endregion LEAF_PARTICIPANT_MAP

#pragma region TREE_GENERATION
    /*
    TREE STRUCTURE:

    |__ Ideas
        |__ BodyEnclosure
        |       |__ Form
        |       |       |__ Rigid1
        |       |       |       |__ Prismatic
        |       |       |       |       |__ Rectangular
        |       |       |       |       |       |__ Multisurface
        |       |       |       |       |       |__ CurvedSurfaces
        |       |       |       |       |       |__ FlatSurfaces
        |       |       |       |       |__ LShaped
        |       |       |       |       |__ Triangular
        |       |       |       |       |__ Hexagon
        |       |       |       |__ Conical
        |       |       |       |__ SemiSpherical
        |       |       |       |__ Cylindrical
        |       |       |               |__ Hourglass
        |       |       |               |__ Assymetrical
        |       |       |               |__ Straight
        |       |       |               |__ SemiCylindrical
        |       |       |       |__ PurseLike
        |       |       |       |__ Sculpted
        |       |       |               |__ ForOtherFunctions
        |       |       |               |__ ToMountAStrap
        |       |       |               |__ ToGraspBody
        |       |       |               |__ ToFormFeet
        |       |       |__ Flexible1
        |       |               |__ Wrappable
        |       |               |__ BeanBag
        |       |__ Dimensions
        |       |__ Colors
        |       |__ Orientation
        |               |__ SelfStanding
        |               |__ VerticalOrientation
        |               |__ HorizontalOrientation
        |               |__ Hangable
    */

    std::shared_ptr<Tree> genealogyTree = std::make_shared<Tree>(Ideas);
    {
        genealogyTree->addNode(BodyEnclosure, Ideas);

        genealogyTree->addNode(Form, BodyEnclosure);
        genealogyTree->addNode(Dimensions, BodyEnclosure);
        genealogyTree->addNode(Colors, BodyEnclosure);
        genealogyTree->addNode(Orientation, BodyEnclosure);

        genealogyTree->addNode(Rigid1, Form);
        genealogyTree->addNode(Flexible1, Form);
        genealogyTree->addNode(SelfStanding, Orientation);
        genealogyTree->addNode(VerticalOrientation, Orientation);
        genealogyTree->addNode(HorizontalOrientation, Orientation);
        genealogyTree->addNode(Hangable, Orientation);

        genealogyTree->addNode(Prismatic, Rigid1);
        genealogyTree->addNode(Conical, Rigid1);
        genealogyTree->addNode(SemiSpherical, Rigid1);
        genealogyTree->addNode(Cylindrical, Rigid1);
        genealogyTree->addNode(PurseLike, Rigid1);
        genealogyTree->addNode(Sculpted, Rigid1);
        genealogyTree->addNode(Wrappable, Flexible1);
        genealogyTree->addNode(BeanBag, Flexible1);

        genealogyTree->addNode(Rectangular, Prismatic);
        genealogyTree->addNode(LShaped, Prismatic);
        genealogyTree->addNode(Triangular, Prismatic);
        genealogyTree->addNode(Hexagon, Prismatic);
        genealogyTree->addNode(Hourglass, Cylindrical);
        genealogyTree->addNode(Assymetrical, Cylindrical);
        genealogyTree->addNode(Straight, Cylindrical);
        genealogyTree->addNode(SemiCylindrical, Cylindrical);
        genealogyTree->addNode(ForOtherFunctions, Sculpted);
        genealogyTree->addNode(ToMountAStrap, Sculpted);
        genealogyTree->addNode(ToGraspBody, Sculpted);
        genealogyTree->addNode(ToFormFeet, Sculpted);

        genealogyTree->addNode(Multisurface, Rectangular);
        genealogyTree->addNode(CurvedSurfaces, Rectangular);
        genealogyTree->addNode(FlatSurfaces, Rectangular);
    }
#pragma endregion TREE_GENERATION

    std::cout << "= INCUBATION TREE =" << std::endl;
    genealogyTree->print();
    std::cout << "-o-o-o-o-o-o-o-o-o-o-o-o-" << std::endl;

    printNoveltyScores("INCUBATION", genealogyTree, leafParticipantMap);
    printVarietyScores("INCUBATION", genealogyTree, levelWeights, leafParticipantMap);
    std::cout << "INCUBATION Number of leaf: " << leafParticipantMap.size() << std::endl;
}

int main()
{
    std::cout << "-TREE SCORE CALCULATOR-" << std::endl;
    CalculateGenealogyTreeScore();
    return 0;
}
