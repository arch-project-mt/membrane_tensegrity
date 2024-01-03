def reformat_obj(obj_file, out_file):
    with open(obj_file, 'r') as f:
        lines = f.readlines()

    with open(out_file, 'w') as f:
        for line in lines:
            if line.startswith('v '):
                f.write(line)
            elif line.startswith('f '):
                new_line = [ls.split('/')[0] for ls in line.split(' ')]
                f.write(" ".join(new_line) + '\n')

if __name__ == '__main__':
    reformat_obj('src/iofiles/donuts.obj', 'src/iofiles/donuts1.obj')
