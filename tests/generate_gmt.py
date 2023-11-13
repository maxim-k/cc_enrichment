import random
import string


def generate_data_string() -> str:
    """
    Generate a data string with:
    - Every line as a term name in format f"Term_{number_of_line}", tab,
      a term description in format f"Description_of_Term_{number_of_line}", tab,
      and a tab separated set of three letter words with the size of the set ranging from 5 to 200.
    - The data has the number of lines ranging from 500 to 2000.
    """

    def random_three_letter_word():
        """Generate a random three-letter word."""
        return "".join(random.choice(string.ascii_uppercase) for _ in range(3))

    # Number of lines to generate
    num_lines = random.randint(500, 2000)

    lines = []
    for i in range(1, num_lines + 1):
        term_name = f"Term_{i}"
        term_description = f"Description_of_Term_{i}"
        num_genes = random.randint(5, 200)
        genes = set()
        while len(genes) < num_genes:
            genes.add(random_three_letter_word())
        line_parts = [term_name, term_description] + list(genes)
        line = "\t".join(line_parts)
        lines.append(line)

    return "\n".join(lines)


# Generating the data string
data_string = generate_data_string()

# Checking the first few lines of the generated data string
data_string.split("\n")[:5]
