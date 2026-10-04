# Build the bundled `forum_posts` dataset: a simulated forum whose sequence
# of topics has second-order memory, used by the hg_sequences()
# documentation and tests.
# Run from the package root: Rscript data-raw/forum_posts.R
#
# Thirty students each write twelve posts. A post is about one of three
# topics (cooking, astronomy, gardening) and holds seven words drawn with
# replacement from that topic's vocabulary and one everyday word shared by
# all topics. The first two topics of a student are drawn at random. After
# two different topics, the next topic is the third one with probability
# 0.8 and otherwise drawn at random; after a repeated topic it is drawn at
# random. Seed 11.

set.seed(11L)
n_students <- 30L
n_turns <- 12L
words_per_post <- 7L
move_probability <- 0.8

vocabulary <- list(
  cooking = c("soup", "onions", "carrots", "simmer", "recipe", "salt",
              "oven", "garlic", "bake", "stew"),
  astronomy = c("telescope", "galaxy", "stars", "planet", "orbit", "comet",
                "nebula", "moon", "eclipse", "astronomers"),
  gardening = c("seeds", "soil", "compost", "tomatoes", "prune", "roses",
                "watering", "mulch", "sprout", "greenhouse")
)
everyday <- c("night", "today", "week", "friends")
topic_names <- names(vocabulary)

# A topic path is sequential (each topic depends on the two before it), so
# it is carried by Reduce().
topic_path <- function() {
  start <- sample(topic_names, 2L, replace = TRUE)
  Reduce(\(path, i) {
    third <- setdiff(topic_names, utils::tail(path, 2L))
    step <- if (length(third) == 1L && stats::runif(1L) < move_probability) {
      third
    } else {
      sample(topic_names, 1L)
    }
    c(path, step)
  }, seq_len(n_turns - 2L), init = start)
}

one_student <- function(student) {
  topic <- topic_path()
  text <- vapply(topic, \(t) {
    paste(c(sample(vocabulary[[t]], words_per_post, replace = TRUE),
            sample(everyday, 1L)), collapse = " ")
  }, character(1L))
  data.frame(post = sprintf("%s_t%02d", student, seq_len(n_turns)),
             student = student, turn = seq_len(n_turns), topic = topic,
             text = unname(text), stringsAsFactors = FALSE)
}

students <- sprintf("st%02d", seq_len(n_students))
forum_posts <- do.call(rbind, lapply(students, one_student))
rownames(forum_posts) <- NULL

stopifnot(
  "one row per post" = nrow(forum_posts) == n_students * n_turns,
  "post ids are unique" = !anyDuplicated(forum_posts$post),
  "no missing values" = !anyNA(forum_posts)
)

save(forum_posts, file = file.path("data", "forum_posts.rda"),
     compress = "xz")
