#[derive(Clone, Debug)]
pub enum GuideTree {
    Leaf(usize),
    Merge {
        left: Box<GuideTree>,
        right: Box<GuideTree>,
    },
}

#[derive(Clone, Debug)]
struct Cluster {
    members: Vec<usize>,
    tree: GuideTree,
}

pub fn build_guide_tree(similarities: &[f32], count: usize) -> GuideTree {
    assert!(count > 0, "guide tree requires at least one sequence");
    if count == 1 {
        return GuideTree::Leaf(0);
    }

    let mut clusters = (0..count)
        .map(|index| Cluster {
            members: vec![index],
            tree: GuideTree::Leaf(index),
        })
        .collect::<Vec<_>>();

    while clusters.len() > 1 {
        let mut best_pair = (0usize, 1usize);
        let mut best_score = f32::NEG_INFINITY;

        for left in 0..clusters.len() {
            for right in (left + 1)..clusters.len() {
                let score = average_similarity(similarities, count, &clusters[left], &clusters[right]);
                if score > best_score {
                    best_score = score;
                    best_pair = (left, right);
                }
            }
        }

        let (left_index, right_index) = best_pair;
        let right = clusters.remove(right_index);
        let left = clusters.remove(left_index);
        let mut merged_members = left.members.clone();
        merged_members.extend_from_slice(&right.members);
        clusters.push(Cluster {
            members: merged_members,
            tree: GuideTree::Merge {
                left: Box::new(left.tree),
                right: Box::new(right.tree),
            },
        });
    }

    clusters.pop().expect("at least one cluster").tree
}

fn average_similarity(similarities: &[f32], sequence_count: usize, left: &Cluster, right: &Cluster) -> f32 {
    let mut total = 0.0f32;
    let mut pairs = 0usize;

    for &i in &left.members {
        for &j in &right.members {
            total += similarities[i * sequence_count + j];
            pairs += 1;
        }
    }

    total / pairs as f32
}
