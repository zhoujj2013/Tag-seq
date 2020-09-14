import svgwrite
import sys
import os

boxWidth = 10
box_size = 15
v_spacing = 3

colors = {'G': '#F5F500', 'A': '#FF5454', 'T': '#00D118', 'C': '#26A8FF', 'N': '#B3B3B3'}

#chr3    46399671        46399698        chr3:46399625-46399726;chr3:46399626-46399727   720     TTCAGGATTCCCGAGTAGCAGATGACC     TTTAGGATTCCCGAGTAGCAGATGACC     CCR5AS;CCR5LB        720;3548
#chr3    46415022        46415049        chr3:46414986-46415087;chr3:46414979-46415080,chr3:46414996-46415097    2389    TTTAGGATTCCCGAGTAGCAGATGACC     TTTAGGATTCCCGAGTAGCAGATGACC  CCR5AS;CCR5LB   2389;4875

def parseSitesFile(infile):
    offtargets = []
    with open(infile, 'r') as f:
        f.readline()
        for line in f:
            line_items = line.split('\t')
            offtarget_sequence = line_items[5]
            offtarget_reads = line_items[4]
            ref_seq = line_items[6]

            names = line_items[7].strip("\n").split(";")
            counts = line_items[8].strip("\n").split(";")

            if offtarget_sequence != '':
                offtargets.append({'seq': offtarget_sequence.strip(),
                                   'reads': int(offtarget_reads.strip()),
                                   'names': names,
                                   'counts': counts
                                   })

    offtargets = sorted(offtargets, key=lambda x: x['reads'], reverse=True)
    return offtargets, ref_seq, names


def visualizeOfftargets(infile, outfile, title):

    output_folder = os.path.dirname("./" + outfile)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get offtargets array from file
    offtargets, ref_seq, names = parseSitesFile(infile)

    # Initiate canvas
    dwg = svgwrite.Drawing(outfile + '.svg', profile='full')

    if title is not None:
        # Define top and left margins
        x_offset = 20
        y_offset = 50
        dwg.add(dwg.text(title, insert=(x_offset, 30), style="font-size:20px; font-family:Courier"))
    else:
        # Define top and left margins
        x_offset = 20
        y_offset = 20

    # Draw ticks
    tick_locations = [1, len(ref_seq)]
    tick_locations += range(len(ref_seq) + 1)[::10][1:]
    for x in tick_locations:
        dwg.add(dwg.text(str(x), insert=(x_offset + (x - 1) * box_size + 2, y_offset - 2), style="font-size:10px; font-family:Courier"))

    # Draw reference sequence row
    for i, c in enumerate(ref_seq):
        y = y_offset
        x = x_offset + i * box_size
        dwg.add(dwg.rect((x, y), (box_size, box_size), fill=colors[c]))
        dwg.add(dwg.text(c, insert=(x + 3, y + box_size - 3), fill='black', style="font-size:15px; font-family:Courier"))

    #dwg.add(dwg.text('Reads', insert=(x_offset + box_size * len(ref_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
    names = title.split(",")
    offset_steps = 0
    for n in names:
        dwg.add(dwg.text(n, insert=(x_offset + box_size * len(ref_seq) + 16 + offset_steps, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))
        offset_steps = offset_steps + 96
    #dwg.add(dwg.text(fm_name_str, insert=(x_offset + box_size * len(ref_seq) + 16, y_offset + box_size - 3), style="font-size:15px; font-family:Courier"))

    # Draw aligned sequence rows
    y_offset += 10  # leave some extra space after the reference row
    copy_x = x_offset
    copy_y = y_offset
    for j, seq in enumerate(offtargets):
        y = y_offset + j * box_size
        copy_y = y_offset + j * box_size
        for i, c in enumerate(seq['seq']):
            x = x_offset + i * box_size
            if c == ref_seq[i] or ref_seq[i] == 'N':
                dwg.add(dwg.text(u"\u2022", insert=(x + 4.5, 2 * box_size + y - 4), fill='black', style="font-size:10px; font-family:Courier"))
            else:
                dwg.add(dwg.rect((x, box_size + y), (box_size, box_size), fill=colors[c]))
                dwg.add(dwg.text(c, insert=(x + 3, 2 * box_size + y - 3), fill='black', style="font-size:15px; font-family:Courier"))
        # draw read counts
        #reads_text = dwg.text(str(seq['reads']), insert=(box_size * (len(ref_seq) + 1) + 20, y_offset + box_size * (j + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
        offset_steps = 0
        for c in seq['counts']:
            reads_text = dwg.text(c, insert=(box_size * (len(ref_seq) + 1) + 20 + offset_steps, y_offset + box_size * (j + 2) - 2), fill='black', style="font-size:15px; font-family:Courier")
            offset_steps = offset_steps + 96
            dwg.add(reads_text)
    
    copy_right_text = dwg.text("Author: zhoujj2013@gmail.com ()", insert=(box_size * (len(ref_seq) + 1) + 20 + offset_steps, copy_y + (offset_steps/96)*15 + 15), fill='black', style="font-size:15px; font-family:Courier")
    dwg.add(copy_right_text)
    blank_line = dwg.text("Dr. Jiajian ZHOU @SMU.China ", insert=(box_size * (len(ref_seq) + 1) + 20 + offset_steps, copy_y + (offset_steps/96)*15 + 15*3), fill='black', style="font-size:15px; font-family:Courier")
    dwg.add(blank_line)

    dwg.save()


def main():
    if len(sys.argv) >= 3:
        if len(sys.argv) == 4:
            title = sys.argv[3]
        else:
            title = sys.argv[2]
        visualizeOfftargets(sys.argv[1], sys.argv[2], title=title)
    else:
        print 'Usage: python visualization.py INFILE out.PREFIX TITLE'


if __name__ == '__main__':
    main()
