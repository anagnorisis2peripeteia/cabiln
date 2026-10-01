// Residue tabs, atom ownership and interaction with one accepted drawing.
// The application supplies selection and highlight policies; this view owns its
// maps and DOM. Clearing a drawing cannot leave another panel's ownership behind.
function createResidueView({ inner, chips: resChips, canHighlight, onSelect }) {
  let residueMap = {}, atomToRes = {}, residueList = [];
  const RES_COLORS = [
    '#2a5080','#2a8050','#802a50','#806a2a','#502a80',
    '#2a6080','#80502a','#2a8070','#6a2a80','#80802a',
    '#3a6080','#3a8060','#603a50','#706a3a','#403a70',
  ];

  let chainData = [];
  let currentBranchSet = new Set();
  let xlinkByRes = {};

  function highlightGroup(idxList) {
    if (!canHighlight()) return;
    clearHighlight();
    activeRIdx = -999;
    const svg = inner.querySelector('svg');
    if (!svg) return;
    svg.classList.add('has-highlight');
    for (const rIdx of idxList) {
      const atoms = residueMap[rIdx] || [];
      for (const aidx of atoms) {
        svg.querySelectorAll(`.atom-${aidx}`).forEach(el =>
          el.classList.add('res-hl'));
      }
    }
    resChips.classList.add('dimmed');
    resChips.querySelectorAll('.res-chip').forEach(c => {
      const ri = parseInt(c.dataset.residue);
      c.classList.toggle('hover', idxList.includes(ri));
    });
    resChips.querySelectorAll('.branch-chip').forEach(c => {
      const cm = JSON.parse(c.dataset.members || '[]');
      // Highlight a branch chip only when ALL its members are in the hover set,
      // not merely "any overlap". Otherwise hovering on one !2 lights up the
      // sibling [!1] and [!3] brackets too because they share the scaffold
      // monomer (e.g. TBMB) — overzealous and confusing for tri-arm scaffolds.
      c.classList.toggle('hover', cm.length > 0 && cm.every(m => idxList.includes(m)));
    });
  }

  function render(data = {}, quality = null) {
    clearHighlight();
    const { residue_map: resMap, residues, layout, crosslink_groups: crosslinkGroups } = data;
    residueMap = resMap || {};
    residueList = residues || [];
    const segments = layout?.segments || [];
    chainData = segments.map(segment => ({ residues: segment.roots }));
    atomToRes = {};
    for (const [rIdx, atoms] of Object.entries(residueMap)) {
      for (const aidx of atoms) atomToRes[aidx] = parseInt(rIdx);
    }

    resChips.innerHTML = '';
    resChips.style.position = '';
    resChips.style.paddingLeft = '';
    currentBranchSet = new Set();
    xlinkByRes = {};
    if (!residueList.length) return;

    const resById = {};
    residueList.forEach((r, i) => { resById[r.idx] = { ...r, colorIdx: i }; });

    const xlinkByMember = {};
    (crosslinkGroups || []).forEach(g => {
      g.members.forEach(mIdx => {
        if (!xlinkByMember[mIdx]) xlinkByMember[mIdx] = [];
        xlinkByMember[mIdx].push(g);
      });
    });
    xlinkByRes = xlinkByMember;

    function makeChip(rIdx, simpleHover) {
      const r = resById[rIdx];
      if (!r) return null;
      const chip = document.createElement('button');
      chip.type = 'button';
      chip.className = 'res-chip';
      chip.textContent = r.abbr;
      chip.dataset.residue = r.idx;
      chip.title = `Select residue ${r.idx + 1}: ${r.abbr}`;
      const assignment = quality?.assignments?.find(item => item.residue_index === r.idx);
      const kind = r.kind || r.quality_kind;
      const preserved = ['opaque', 'synthetic'].includes(kind) ? kind :
        !kind && (r.is_synthetic || assignment?.recognized === false) ? 'synthetic' : '';
      if (preserved) {
        chip.dataset.quality = preserved;
        chip.title += ` · ${preserved === 'opaque' ? 'Opaque preserved fragment' : 'Synthetic preserved region'}; available sites remain editable`;
        chip.setAttribute('aria-label', chip.title);
      }
      chip.style.background = RES_COLORS[r.colorIdx % RES_COLORS.length];
      const xlinks = xlinkByMember[r.idx];
      if (!simpleHover && xlinks && xlinks.length) {
        const allMembers = [...new Set(xlinks.flatMap(g => g.members))];
        chip.addEventListener('mouseenter', () => view.highlightGroup(allMembers));
      } else {
        chip.addEventListener('mouseenter', () => view.highlightResidue(r.idx));
      }
      chip.addEventListener('mouseleave', clearHighlight);
      chip.addEventListener('focus', () => view.highlightResidue(r.idx));
      chip.addEventListener('blur', clearHighlight);
      chip.addEventListener('click', () => onSelect(r));
      return chip;
    }

    function makeSeparator(text, memberIdxs) {
      const el = document.createElement('span');
      el.className = 'res-chip branch-chip';
      el.style.background = '#3a3a50';
      el.style.fontWeight = '700';
      el.textContent = text;
      el.dataset.members = JSON.stringify(memberIdxs || []);
      if (memberIdxs && memberIdxs.length) {
        el.tabIndex = 0;
        el.setAttribute('aria-label', `Highlight ${text} group`);
        el.addEventListener('mouseenter', () => view.highlightGroup(memberIdxs));
        el.addEventListener('mouseleave', clearHighlight);
        el.addEventListener('focus', () => view.highlightGroup(memberIdxs));
        el.addEventListener('blur', clearHighlight);
      }
      return el;
    }

    function makeXlinkChip(tag, members) {
      const el = document.createElement('span');
      el.className = 'res-chip branch-chip xlink-chip';
      el.style.background = '#503a4a';
      el.style.fontWeight = '700';
      el.style.fontSize = '0.8em';
      el.textContent = tag;
      el.dataset.members = JSON.stringify(members);
      el.tabIndex = 0;
      el.setAttribute('aria-label', `Highlight connection ${tag}`);
      el.addEventListener('mouseenter', () => view.highlightGroup(members));
      el.addEventListener('mouseleave', clearHighlight);
      el.addEventListener('focus', () => view.highlightGroup(members));
      el.addEventListener('blur', clearHighlight);
      return el;
    }

    const groups = layout?.groups || [];
    const markers = layout?.markers || [];
    function indexBy(items, field) {
      const index = new Map();
      for (const item of items) {
        if (!index.has(item[field])) index.set(item[field], []);
        index.get(item[field]).push(item);
      }
      return index;
    }
    const groupsByHost = indexBy(groups, 'host');
    const markersByGroup = indexBy(markers, 'group');
    const markersByResidue = indexBy(markers, 'residue');
    currentBranchSet = new Set(groups.flatMap(group => group.members));

    function expandWithXlinks(ids) {
      const members = new Set(ids);
      for (const id of ids) {
        for (const link of xlinkByMember[id] || []) {
          link.members.forEach(member => members.add(member));
        }
      }
      return [...members];
    }

    function groupMembers(group) {
      return expandWithXlinks([
        ...group.members,
        ...(markersByGroup.get(group.id) || []).flatMap(marker => marker.members),
      ]);
    }

    function appendGroup(group) {
      const members = groupMembers(group);
      if (group.opening) resChips.appendChild(makeSeparator(group.opening, members));
      // A marker-only bracket belongs to its host but contains no monomer.
      if (!group.roots.length) {
        (markersByGroup.get(group.id) || []).forEach(marker => {
          resChips.appendChild(makeXlinkChip(marker.tag, marker.members));
        });
      }
      group.roots.forEach(id => appendOccurrence(id, group.id));
      if (group.closing) resChips.appendChild(makeSeparator(group.closing, members));
    }

    function appendOccurrence(id, context = null) {
      const ownMarkers = (markersByResidue.get(id) || []).filter(marker => marker.group === context);
      ownMarkers.filter(marker => marker.before).forEach(marker => {
        resChips.appendChild(makeXlinkChip(marker.tag, marker.members));
      });
      const links = xlinkByMember[id] || [];
      const chip = makeChip(id, context !== null || links.length > 1);
      if (chip) resChips.appendChild(chip);
      const children = (groupsByHost.get(id) || []).filter(group => group.parent === context);
      if (children.length || links.length > 1) {
        const members = expandWithXlinks([id, ...children.flatMap(groupMembers)]);
        const whole = makeSeparator('$', members);
        whole.style.background = '#3a5050';
        resChips.appendChild(whole);
      }
      children.forEach(appendGroup);
      ownMarkers.filter(marker => !marker.before).forEach(marker => {
        resChips.appendChild(makeXlinkChip(marker.tag, marker.members));
      });
    }

    for (const segment of segments) {
      if (segments.length > 1) {
        resChips.appendChild(makeSeparator('%', expandWithXlinks(segment.members)));
      }
      segment.roots.forEach(id => appendOccurrence(id));
    }

    wireUpSvgHover();
  }

  let activeRIdx = null;

  function highlightResidue(rIdx) {
    if (!canHighlight()) return;
    if (rIdx === activeRIdx) return;
    clearHighlight();
    activeRIdx = rIdx;
    const atoms = residueMap[rIdx] || [];
    const svg = inner.querySelector('svg');
    if (!svg) return;

    svg.classList.add('has-highlight');
    for (const aidx of atoms) {
      svg.querySelectorAll(`.atom-${aidx}`).forEach(el =>
        el.classList.add('res-hl'));
    }

    resChips.classList.add('dimmed');
    resChips.querySelectorAll('.res-chip').forEach(c =>
      c.classList.toggle('hover', parseInt(c.dataset.residue) === rIdx));
  }

  function clearHighlight() {
    activeRIdx = null;
    const svg = inner.querySelector('svg');
    if (svg) {
      svg.classList.remove('has-highlight');
      svg.querySelectorAll('.res-hl').forEach(el => el.classList.remove('res-hl'));
    }
    resChips.classList.remove('dimmed');
    resChips.querySelectorAll('.res-chip.hover').forEach(c => c.classList.remove('hover'));
  }

  function svgResidueIndex(target, svg) {
    for (let el = target; el && el !== svg; el = el.parentElement) {
      const atom = (el.getAttribute('class') || '').match(/atom-(\d+)/);
      if (atom) {
        const rIdx = atomToRes[parseInt(atom[1])];
        if (rIdx !== undefined) return rIdx;
      }
    }
  }

  function wireUpSvgHover() {
    const svg = inner.querySelector('svg');
    if (!svg) return;
    svg.addEventListener('mousemove', e => {
      if (!canHighlight()) return;
      const rIdx = svgResidueIndex(e.target, svg);
      if (rIdx === undefined) return clearHighlight();
      const xlinks = xlinkByRes[rIdx];
      if (xlinks && xlinks.length) {
        view.highlightGroup([...new Set(xlinks.flatMap(g => g.members))]);
      } else {
        view.highlightResidue(rIdx);
      }
    });
    svg.addEventListener('mouseleave', clearHighlight);
    svg.addEventListener('click', e => {
      const rIdx = svgResidueIndex(e.target, svg);
      const residue = residueList.find(r => r.idx === rIdx);
      if (residue) onSelect(residue);
    });
  }


  function clear() { render(); }

  function setStale(stale) {
    resChips.classList.toggle('stale', stale);
    resChips.setAttribute('aria-disabled', String(stale));
    resChips.querySelectorAll('button').forEach(button => { button.disabled = stale; });
    if (stale) clearHighlight();
  }

  function select(left = null, right = null) {
    resChips.querySelectorAll('.res-chip').forEach(chip => {
      const index = parseInt(chip.dataset.residue);
      chip.style.outline = index === left ? '2px solid #5a9ae0'
        : index === right ? '2px solid #e0a05a' : '';
    });
  }

  function insertionAnchor(left, right) {
    if (left == null || right == null || left === right ||
        currentBranchSet.has(left) || currentBranchSet.has(right)) return null;
    const chain = chainData.find(item => item.residues.includes(left) && item.residues.includes(right));
    if (!chain) return null;
    const a = chain.residues.indexOf(left), b = chain.residues.indexOf(right);
    return Math.abs(a - b) === 1 ? chain.residues[Math.min(a, b)] : null;
  }

  const view = {
    render, clear, setStale, select, insertionAnchor,
    clearHighlight, highlightGroup, highlightResidue,
    residue: index => residueList.find(residue => residue.idx === index),
    get residues() { return residueList; },
    get atoms() { return residueMap; },
    get atomOwners() { return atomToRes; },
  };
  return view;
}
