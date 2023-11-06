import sys

def count_question_marks(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
    
    question_mark_count = content.count('?')
    
    return question_mark_count

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python count_question_marks.py file1.txt file2.txt ...")
        sys.exit(1)
    
    file_paths = sys.argv[1:]
    
    total_question_mark_count = 0
    
    for file_path in file_paths:
        question_mark_count = count_question_marks(file_path)
        print(f"File: {file_path} - Question Marks: {question_mark_count}")
        total_question_mark_count += question_mark_count
    
    print(f"Total Question Marks in all Files: {total_question_mark_count}")
